import subprocess as sp
import matplotlib.pyplot as plt
import cc3d
import numpy as np
import sys
from copy import deepcopy
from pathlib import Path

import trimesh
import pymeshfix
import ipywidgets as widgets

import matplotlib.image as mpimg
import nibabel as nib
from nibabel.freesurfer.io import read_geometry, write_geometry
from skimage import measure
from parse import parse
from tqdm.notebook import tqdm
from PIL import Image, ImageFont, ImageDraw
import pandas as pd
import signal

import sys
import csv
import shutil

import mne
from mne.surface import decimate_surface, _triangle_neighbors
from mne.bem import _get_solids
from mne.io.constants import FIFF
from mne.viz import plot_alignment, set_3d_view

from eegip.atlas import Atlas
from eegip import DatasetMng

#out = widgets.Output(layout={'border': '1px solid black'})


def run_command(command, dry_run=False):
    print("------------------ RUNNING A COMMAND IN BASH ------------------")
    print(command)    
    print("-------------------------- OUTPUT -----------------------------")
    if not dry_run:
        bash(command.replace("\\", "").split())
    print("------------------------ END OUTPUT ---------------------------\n")    


class VerboseCalledProcessError(sp.CalledProcessError):
    def __str__(self):
        if self.returncode and self.returncode < 0:
            try:
                msg = "Command '%s' died with %r." % (
                    self.cmd, signal.Signals(-self.returncode))
            except ValueError:
                msg = "Command '%s' died with unknown signal %d." % (
                    self.cmd, -self.returncode)
        else:
            msg = "Command '%s' returned non-zero exit status %d." % (
                self.cmd, self.returncode)

        return f'{msg}\n' \
               f'Stdout:\n' \
               f'{self.output}\n' \
               f'Stderr:\n' \
               f'{self.stderr}'


def bash(cmd, print_stdout=True, print_stderr=True):
    proc = sp.Popen(cmd, stderr=sp.PIPE, stdout=sp.PIPE)

    all_stdout = []
    all_stderr = []
    while proc.poll() is None:
        for stdout_line in proc.stdout:
            if stdout_line != '':
                if print_stdout:
                    print(stdout_line.decode(), end='')
                all_stdout.append(stdout_line)
        for stderr_line in proc.stderr:
            if stderr_line != '':
                if print_stderr:
                    print(stderr_line.decode(), end='', file=sys.stderr)
                all_stderr.append(stderr_line)

    stdout_text = ''.join([x.decode() for x in all_stdout])
    stderr_text = ''.join([x.decode() for x in all_stderr])
    if proc.wait() != 0:
        raise VerboseCalledProcessError(proc.returncode, cmd, stdout_text, stderr_text)


def check_bem(fs_subject_path, labels_vol, subject_):
    atlas_ = Atlas(subjects_dir=fs_subject_path, fs_subject=subject_)
    atlas_.build_bem(labels_vol=labels_vol)
    atlas_.mesh_atlas_pacels()
    atlas_.compute_parcels_centers_of_masse()


# In some of our MRI, there appears to have artifacts that shows like small line segments
# over the background. We correct these by zeroing any small separated clusted of non null 
# voxels.    
def correct_line_artefact(epi_img_data_):
    components = cc3d.connected_components((epi_img_data_ != 0).astype(int))
    label_id, count = np.unique(components, return_counts=True)
    id_zero, id_non_zero = label_id[count > 500000]
    ind_artefact = np.where((components != id_zero) & (components != id_non_zero))
    epi_img_data_[ind_artefact] = 0

    
def process_bem(fs_subject_path, bem_path_):
    """
     Main function for BEM extraction.
    """
    subject_ = bem_path_.name.split("_")[0].replace("AVG", "ANTS")

    epi_img_ = nib.load(str(bem_path_))

    for file_name_, discard_inds in zip(["outer_skin.surf", "outer_skull.surf", "inner_skull.surf"],
                                        [[1, 2, 3, 4], [1, 2, 3], [1, 2]]):
        epi_img_data_ = deepcopy(epi_img_.get_fdata())

        correct_line_artefact(epi_img_data_)

        cond = np.stack([(epi_img_data_ == i_) for i_ in discard_inds]).sum(0)
        epi_img_data_[np.where(cond)] = 1
        epi_img_data_[np.where(np.logical_not(cond))] = 0
        vertices, simplices = measure.marching_cubes_lewiner(epi_img_data_, spacing=(1, 1, 1),
                                                             allow_degenerate=False)[:2]

        path_white = fs_subject_path / subject_ / "surf" / "lh.white"

        try:
            volume_info_ = read_geometry(path_white, read_metadata=True)[2]
        except:
            print("Skipping subject {}...".format(subject_))
            continue

        vertices = vertices @ epi_img_.affine[:3, :3] + epi_img_.affine[:3, 3] - volume_info_["cras"]

        mesh_ = trimesh.Trimesh(vertices=vertices, faces=simplices)
        trimesh.repair.fix_normals(mesh_, multibody=False)

        smooth_mesh = trimesh.smoothing.filter_laplacian(deepcopy(mesh_), lamb=0.8, iterations=15,
                                                         volume_constraint=True)

        bem_output_path_ = fs_subject_path / subject_ / "bem"
        bem_output_path_.mkdir(parents=True, exist_ok=True)

        vertices, faces_ = smooth_mesh.vertices, smooth_mesh.faces

        # Defect corrections for the large meshes
        vertices, faces_ = fix_all_defects(vertices, faces_)

        # Writing a freesufer mesh file
        file_name_large = file_name_.split(".")[0] + "_large.surf"
        write_geometry(str(bem_output_path_ / file_name_large),
                       vertices, faces_)

        # Writing an obj mesh file
        with (bem_output_path_ / file_name_large).with_suffix(".obj").open('w') as file_obj:
            file_obj.write(trimesh.exchange.obj.export_obj(trimesh.Trimesh(vertices, faces_)))

        # Decimating BEM surfaces
        vertices, faces_ = decimate_surface(vertices, faces_,
                                            n_triangles=5120, method='sphere')

        # Defect correction for decimated meshes...
        vertices, faces_ = fix_all_defects(vertices, faces_)

        # Writing an obj mesh file
        with (bem_output_path_ / file_name_).with_suffix(".obj").open('w') as file_obj:
            file_obj.write(trimesh.exchange.obj.export_obj(trimesh.Trimesh(vertices, faces_)))

        # Writing a freesufer mesh file
        print("Writing {}...".format(str(bem_output_path_ / file_name_)))
        write_geometry(str(bem_output_path_ / file_name_),
                       vertices, faces_, volume_info=volume_info_)
        

def check_mesh(vertices, faces_):
    assert (surface_is_complete(vertices, faces_))
    assert (not has_topological_defects(vertices, faces_))
    assert (not has_degenerated_faces(vertices, faces_))
    assert trimesh.Trimesh(vertices, faces_).is_watertight     
        
        
def fix_all_defects(vertices, faces_):   
    if has_degenerated_faces(vertices, faces_):
        vertices, faces_ = remove_degenerated_faces(vertices, faces_)
        assert (not has_degenerated_faces(vertices, faces_))

    if has_topological_defects(vertices, faces_):
        print("The decimated mesh has topological defects. Fixing it.")
        vertices, faces_ = fix_topological_defects(vertices, faces_)
        if has_degenerated_faces(vertices, faces_):
            vertices, faces_ = remove_degenerated_faces(vertices, faces_)
        assert (not has_topological_defects(vertices, faces_))

    if not surface_is_complete(vertices, faces_) or not trimesh.Trimesh(vertices, faces_).is_watertight:
        print("The decimated mesh has holes. Fixing it.")
        vertices, faces_ = repair_holes(vertices, faces_)

    check_mesh(vertices, faces_)

    return vertices, faces_      
            
        
def surface_is_complete(vertices, faces_):
    """Check the sum of solid angles as seen from inside."""
    cm = vertices.mean(axis=0)
    tot_angle = _get_solids(vertices[faces_], cm[np.newaxis, :])[0]
    prop = tot_angle / (2 * np.pi)
    return np.abs(prop - 1.0) < 1e-5


def correction_two_neighboring_tri(vertices, faces_, faulty_vert_ind):
    ind_faces_to_remove = []
    new_faces = []
    for ind in faulty_vert_ind:
        ind_faces = np.where(faces_ == ind)[0]
        ind_faces_to_remove.extend(ind_faces)
        face1, face2 = faces_[ind_faces]
        new_face = np.unique(np.concatenate((face1, face2)))
        new_face = np.delete(new_face, np.where(new_face == ind))
        assert (len(new_face) == 3)  # If == 4, it means that face1 and face2 do not share a common edge

        new_det = np.linalg.det(vertices[new_face])
        assert new_det  # If zero, the three points are colinear

        # Align the normals
        det1 = np.linalg.det(vertices[face1])
        if np.sign(det1) == np.sign(new_det):
            new_face = new_face[[1, 0, 2]]

        new_faces.append(new_face)

    return np.array(ind_faces_to_remove, dtype=int), new_faces
    
    
def reindex_vertices(vertices, faces_, ind_vertices_to_remove):
    decrement = np.cumsum(np.zeros(vertices.shape[0], dtype=int) +
                          np.in1d(np.arange(vertices.shape[0]), ind_vertices_to_remove))
    vertices = np.delete(vertices, ind_vertices_to_remove, axis=0)
    faces_ = faces_ - decrement[faces_]
    return vertices, faces_
    
    
def get_topological_defects(vertices, faces_):
    #    Find neighboring triangles, accumulate vertex normals, normalize
    neighbor_tri = _triangle_neighbors(faces_, len(vertices))

    #   Check for topological defects
    zero, one, two = list(), list(), list()
    for ni, n in enumerate(neighbor_tri):
        if len(n) < 3:
            if len(n) == 0:
                zero.append(ni)
            elif len(n) == 1:
                one.append(ni)
            else:
                two.append(ni)

    return zero, one, two
    
    
def has_topological_defects(vertices, faces_):
    zero, one, two = get_topological_defects(vertices, faces_)
    return len(zero) or len(one) or len(two)
        
    
# Code extracted and slighly modified from mne.surface.complete_surface_info
# for compactness of the example
def fix_topological_defects(vertices, faces_):
    zero, one, two = get_topological_defects(vertices, faces_)
    ind_faces_to_remove = []
    if len(zero) > 0:
        print('    Vertices do not have any neighboring '
              'triangles: [%s]' % ', '.join(str(z) for z in zero))
        print('    Correcting by removing these vertices.')

    if len(one) > 0:
        print('    Vertices have only one neighboring '
              'triangles: [%s]'
              % ', '.join(str(tri) for tri in one))
        print('    Correcting by removing these vertices and their neighboring triangles.')
        ind_faces_to_remove.extend(np.where(faces_ == one)[0].tolist())

    if len(two) > 0:
        print('    Vertices have only two neighboring '
              'triangles, removing neighbors: [%s]'
              % ', '.join(str(tri) for tri in two))
        print('    Correcting by merging the two neighboring '
              'triangles and removing the faulty vertices.')

        ind_faces, faces_to_add = correction_two_neighboring_tri(vertices, faces_, two)
        ind_faces_to_remove.extend(ind_faces)
        faces_ = np.concatenate((np.delete(faces_, np.array(ind_faces_to_remove, dtype=int), axis=0), faces_to_add))

    vertices_to_remove = np.concatenate((zero, one, two)).astype(int)
    if len(vertices_to_remove):
        vertices, faces_ = reindex_vertices(vertices, faces_, vertices_to_remove)
    else:
        print("No issue found with the mesh.")

    return vertices, faces_
        

def has_degenerated_faces(vertices, faces_):
    return not np.all(trimesh.Trimesh(vertices, faces_).remove_degenerate_faces())


def remove_degenerated_faces(vertices, faces_):
    mesh_ = trimesh.Trimesh(vertices, faces_)
    mesh_.remove_degenerate_faces()
    return mesh_.vertices, mesh_.faces  

    
def repair_holes(vertices, faces_):
    # trimesh has a hole fixing function, but it just deals with
    # 3 or 4 vertices holes.
    meshfix = pymeshfix.MeshFix(vertices, faces_)
    meshfix.repair()
    vertices = meshfix.v  # numpy np.float array
    faces_ = meshfix.f  # numpy np.int32 array

    # The return mesh has a solid angle of -1 instead of 1.
    # Correcting this.
    mesh_ = trimesh.Trimesh(vertices=vertices, faces=faces_)
    trimesh.repair.fix_normals(mesh_, multibody=False)
    return mesh_.vertices, mesh_.faces


def show_slices(slices):
    """ Function to display row of image slices """
    axes_ = plt.subplots(1, len(slices), figsize=(12, 5))[1]
    for slice_, ax_ in zip(slices, axes_):
        ax_.imshow(slice_.T, cmap="gray", origin="lower")
        

def save_sfp_montage(montage, fpath):
    assert (montage.dig[0]["kind"] == FIFF.FIFFV_POINT_CARDINAL and
            montage.dig[0]["ident"] == FIFF.FIFFV_POINT_LPA)
    assert (montage.dig[1]["kind"] == FIFF.FIFFV_POINT_CARDINAL and
            montage.dig[1]["ident"] == FIFF.FIFFV_POINT_NASION)
    assert (montage.dig[2]["kind"] == FIFF.FIFFV_POINT_CARDINAL and
            montage.dig[2]["ident"] == FIFF.FIFFV_POINT_RPA)

    with fpath.open("w") as file_obj:
        for dig_, name in zip(montage.dig[:3], ["FidT9", "FidNz", "FidT10"]):
            file_obj.write("{}\t{}\t{}\t{}\n".format(name, dig_["r"][0], dig_["r"][1], dig_["r"][2]))
        for dig_, name in zip(montage.dig[3:], montage.ch_names):
            file_obj.write("{}\t{}\t{}\t{}\n".format(name, dig_["r"][0], dig_["r"][1], dig_["r"][2]))


def surface_intersect(inner_surface_path, outer_surface_path):

    inner_mesh = trimesh.Trimesh(*read_geometry(inner_surface_path, read_metadata=False)) 
    outer_mesh = trimesh.Trimesh(*read_geometry(outer_surface_path, read_metadata=False))  

    ray_interceptor_inner = trimesh.ray.ray_pyembree.RayMeshIntersector(inner_mesh)
    inds = ray_interceptor_inner.intersects_location(ray_origins=outer_mesh.vertices, 
                                                     ray_directions=outer_mesh.vertex_normals,
                                                     multiple_hits=True)[1]    
    return len(inds) != 0


def get_source_file(mnd_root_path, subject, space, volume_type, extension="nii.gz", missing="raise"):
    for subset in ["Infants", "Preschool"]:
        file_name = f"{subject}_{volume_type}.{extension}"
        path = mnd_root_path / subset / space / file_name
        
        if path.exists():
            return path        
        
    if missing == "raise":
        raise FileNotFoundError("No source file {} available for subject {}.".format(str(path), subject))
    return path


def import_volumes(subject, fs_subject_path, mnd_root_path):
    command = ("mri_convert --conform {}\\\n" + 
               "                      {}")    
 
    for volume_type, space in zip(["brain", "head", "head_brain"], ["Brain", "Head", "Head"]):
        run_command(command.format(get_source_file(mnd_root_path, subject, space, volume_type),
                                   fs_subject_path / subject / f"{volume_type}.nii.gz"))
        
    for volume_type in ["brain_ANTS_IXI_atlas", "brain_ANTS_LPBA40_atlas", 
                        "brain_IXI_atlas", "brain_LPBA40_atlas"]:
        try:
            run_command(command.format(get_source_file(mnd_root_path, subject, "Brain/Atlas", volume_type),
                                       fs_subject_path / subject / f"{volume_type}.nii.gz"))
        except FileNotFoundError:
            continue       
        
        
def make_brain_mask(subject, mnd_root_path, fs_subject_path):
    command = "mri_binarize --i {} \\\n             --o {}  --min 1 --max 1"
    bem_path = get_source_file(mnd_root_path, subject.replace("ANTS", "AVG"), "Sources/BEM", "segmented_BEM4")
    head_brain_mask_path = fs_subject_path / subject / "head_brain_mask.nii.gz" 
    run_command(command.format(bem_path, head_brain_mask_path))
    
    
def create_flirt_mask(subject, fs_subject_path):
    command = "mri_binarize --i {} \\\n             --o {}  --min 1"
    brain_path = fs_subject_path / subject / "brain.nii.gz"
    brain_mask_path = fs_subject_path / subject / "brain_mask.nii.gz"        
    run_command(command.format(brain_path, brain_mask_path))


def map_brain_to_head_space(subject, fs_subject_path):
    brain_path = fs_subject_path / subject / "brain.nii.gz"
    head_brain_path = fs_subject_path / subject / "head_brain.nii.gz"

    brain_mask_path = fs_subject_path / subject / "brain_mask.nii.gz"

    brain_out_path = fs_subject_path / subject / "mprage.nii.gz"                 
    trans_mat_path = fs_subject_path / subject / "flirt_trans_mat.txt"            

    # By not using a weight mask on the ref surface, we want to force the brain to fit within
    # the brain volume defined in the head space. If it does not, the brain meshes cannot be fit 
    # within the bem surfaces of the head. 
    # This seems particularly important for the yougest template    
    command = ("flirt -dof 12 -setbackground 0 \\\n" +
               "      -in {} \\\n" +
               "      -inweight {} \\\n" +
               "      -ref {} \\\n" +              # "      -refweight {} \\\n" +
               "      -omat {} \\\n" +
               "      -out {}")

    command = command.format(brain_path, brain_mask_path, 
                             head_brain_path,
                             trans_mat_path, brain_out_path)                 
    run_command(command)    
    

def map_atlases(subject, fs_subject_path):
    head_brain_path = fs_subject_path / subject / "head_brain.nii.gz"
    trans_mat_path = fs_subject_path / subject / "flirt_trans_mat.txt"            

    for atlas in ["brain_ANTS_IXI_atlas", "brain_ANTS_LPBA40_atlas", 
                  "brain_IXI_atlas", "brain_LPBA40_atlas"]:

        command = ("flirt -applyxfm \\\n" + 
                   "      -in {} \\\n" +
                   "      -ref {} \\\n" +   
                   "      -init {} \\\n" +
                   "      -out {}")

        command = command.format(fs_subject_path / subject / "{}.nii.gz".format(atlas), 
                                 head_brain_path, trans_mat_path, 
                                 fs_subject_path / subject / "{}_head_space.nii.gz".format(atlas))                 
        run_command(command)    

    
def create_recon_mask(subject, fs_subject_path):
    command = "mri_nu_correct.mni --i {} \\\n                   --o {} --n 2"
    command = command.format(str(fs_subject_path / subject / "mprage.nii.gz"), 
                             str(fs_subject_path / subject / "brainmask_nu.nii.gz"))
    run_command(command)


def adding_T1(subject, fs_subject_path):
    # Using the head volume as the T1 because this is used
    # for things like general visualization and fiducial localization
    # and having the head structure (not just the skull striped brain)
    # is useful for that purpose.
    orig = str(fs_subject_path / subject / "head.nii.gz")
    target = fs_subject_path / subject / "mri" / "T1.mgz"

    command = ("mri_convert {orig} \\\n"
               "            {target}").format(orig=orig, target=str(target))
    target.parent.mkdir(parents=True, exist_ok=True)
    run_command(command)


def recon_all(subject, fs_subject_path):
    age = parse("ANTS{}3T", subject).fixed[0]

    if "Weeks" in age:
        command = ('infant_recon_all --s {subject} --newborn \\\n' +
                   '                 --masked {path}/brainmask_nu.nii.gz')
        command = command.format(subject=subject, path=str(fs_subject_path / subject))
    else:
        if "Months" in age:
            age_in_months = float(age[:-6].replace("-", "."))
        elif "Years" in age:
            age_in_months = float(age[:-5].replace("-", ".")) * 12
        else:
            raise ValueError

        command = ('infant_recon_all --s {subject} --age {age} \\\n' +
                   '                 --masked {path}/brainmask_nu.nii.gz')
        command = command.format(subject=subject, age=int(age_in_months),
                                 path=str(fs_subject_path / subject))

    run_command(command)


def extract_brain_surfaces(subject, fs_subject_path, mnd_root_path):
    # Making the FreeSurfer subject folders
    (fs_subject_path / subject).mkdir(exist_ok=True)

    # Commands for converting the volumes in the "conform" space
    import_volumes(subject, fs_subject_path, mnd_root_path)

    # Making a brain mask for the head space, using the BEM
    make_brain_mask(subject, mnd_root_path, fs_subject_path)

    # Creating flirt masks
    create_flirt_mask(subject, fs_subject_path)

    # FLIRT transform of the brain to the head space
    map_brain_to_head_space(subject, fs_subject_path)

    # FLIRT transform the IXI and LPBA40 atlases
    map_atlases(subject, fs_subject_path)

    # Commands for running intensity normalization for the brain masks for infant-recon-all
    # Since the brain is already skull stripped, we just use the brain as mask.
    create_recon_mask(subject, fs_subject_path)

    # Commands for adding the mri/T1.mgz file expected by MNE
    adding_T1(subject, fs_subject_path)

    # Commands to run infant_recon_all
    recon_all(subject, fs_subject_path)


def fix_intersecting_surfaces(inner_surface_path, outer_surface_path,
                              move_margin=1.5, out_path=None):

    inner_mesh = trimesh.Trimesh(*read_geometry(inner_surface_path, read_metadata=False)) 
    outer_mesh = trimesh.Trimesh(*read_geometry(outer_surface_path, read_metadata=False))  

    edges = np.array(trimesh.graph.vertex_adjacency_graph(outer_mesh).edges)
        
    ray_interceptor_inner = trimesh.ray.ray_pyembree.RayMeshIntersector(inner_mesh)

    outer_surface_path = Path(outer_surface_path)
    stem = outer_surface_path.stem
    name = outer_surface_path.name
    if out_path is None:
        out_path = outer_surface_path.parent / (stem + suffix + name[len(stem):])    
    
    out_vertices = deepcopy(outer_mesh.vertices)    
    
    # Instead of using the vertex normal, we take the median of the 
    # normal of the neighbor vertex in order for this direction
    # to be more robust to sharp angular changes.
    x = np.concatenate([edges[:, [1, 0]], edges])
    normal = [np.median(outer_mesh.vertex_normals[x[x[:, 0] == i][:, 1]], 0) 
              for i in np.arange(outer_mesh.vertices.shape[0])] 
    
    # Intersecting surface correction (from outer to inner)
    intersections, inds = ray_interceptor_inner.intersects_location(ray_origins=outer_mesh.vertices, 
                                                                    ray_directions=normal,
                                                                    multiple_hits=True)[:2]    
    delta1 = np.zeros(out_vertices.shape[0])
    
    print("inner: ", inner_surface_path)
    print("outer: ", outer_surface_path)

    if len(inds):
        dist = np.sqrt(((outer_mesh.vertices[inds, :]-intersections)**2).sum(axis=1))
        # <5 to avoid picking very long distance points        
        inds = inds[dist < 5]   
        dist = dist[dist < 5]
        if len(inds):
            msg = "The inner surface intersect the outer surface. Pushing back the outer " + \
                  "surface {} mm out of the inner surface. Saving the outer ".format(move_margin) + \
                  "surface as {}.".format(out_path)
            print(msg)
            # <5 to avoid picking very long distance points
            delta1[inds] = dist + move_margin

    # Intersecting surface correction (from inner to outer)
    delta2 = np.zeros(out_vertices.shape[0])

    closest, distance, triangle_id = trimesh.proximity.closest_point(outer_mesh, inner_mesh.vertices)

    normal = (inner_mesh.vertices - closest)/distance[:, None]
    face_normals = outer_mesh.face_normals[triangle_id]

    angles = trimesh.geometry.vector_angle(np.stack([normal, face_normals], axis=1))

    new_deltas_df = pd.DataFrame(dict(vertex_id=outer_mesh.faces[triangle_id[angles < 1.0]].ravel(), 
                                      delta=np.tile(distance[angles < 1.0] + move_margin, [3, 1]).T.ravel())
                                 ).groupby("vertex_id").max().reset_index()

    delta2[new_deltas_df.vertex_id] = new_deltas_df.delta    

    deltas = np.stack([delta1, delta2]).max(axis=0)    
    
    if np.all(deltas == 0):
        return
        
    out_vertices += deltas[:, None]*outer_mesh.vertex_normals

    volume_info = read_geometry(outer_surface_path, read_metadata=True)[2]
    write_geometry(out_path, out_vertices, outer_mesh.faces, volume_info=volume_info)
    return deltas


def correct_intersecting_meshes(fs_subject_path, subject, suffix=""):
    deltas = {}
    pial_lh_path = fs_subject_path / subject / "surf" / "lh.pial"
    pial_rh_path = fs_subject_path / subject / "surf" / "rh.pial"

    inner_skull_path_pattern = fs_subject_path / subject / "bem" / 'inner_skull{}.surf'
    outer_skull_path_pattern = fs_subject_path / subject / "bem" / 'outer_skull{}.surf'
    outer_skin_path_pattern = fs_subject_path / subject / "bem" / 'outer_skin{}.surf'       
    if suffix != "":  
        for path_pattern in [inner_skull_path_pattern, outer_skull_path_pattern, outer_skin_path_pattern]:
            shutil.copy(str(path_pattern).format(""), str(path_pattern).format(suffix))

    inner_skull_path = Path(str(inner_skull_path_pattern).format(suffix))
    outer_skull_path = Path(str(outer_skull_path_pattern).format(suffix))
    outer_skin_path = Path(str(outer_skin_path_pattern).format(suffix))
        
    inner_paths = [pial_lh_path, pial_rh_path, inner_skull_path, outer_skull_path]
    outer_paths = [inner_skull_path, inner_skull_path, outer_skull_path, outer_skin_path] 
    for inner_path, outer_path in zip(inner_paths, outer_paths):
        deltas[inner_path.name] = fix_intersecting_surfaces(inner_path, outer_path, out_path=outer_path, 
                                                            move_margin=0.5)
    return deltas
