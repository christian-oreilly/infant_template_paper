import os
import mne
from pathlib import Path
from mayavi import mlab
from mne.viz import plot_alignment, set_3d_view
from tqdm import tqdm
import numpy as np

os.environ["ETS_TOOLKIT"] = "wx"
fs_subject_path = Path("/usr/local/freesurfer/subjects/")

path = Path("/home/christian/synchedin/infants_atlas_modeling/infant_template_paper")


subjects = ["ANTS2-0Weeks3T", "ANTS1-0Months3T", "ANTS2-0Months3T", "ANTS3-0Months3T", 
            "ANTS4-5Months3T", "ANTS6-0Months3T", "ANTS7-5Months3T", "ANTS9-0Months3T", 
            "ANTS10-5Months3T", "ANTS12-0Months3T", "ANTS15-0Months3T", "ANTS18-0Months3T", 
            "ANTS2-0Years3T"]

def save_images(subject_, montage_fname_):
    montage = mne.channels.read_dig_fif(montage_fname_)    
    info = mne.create_info(ch_names=montage.ch_names, ch_types=['eeg']*len(montage.ch_names), sfreq=100.0)
    info.set_montage(montage)
    trans = mne.channels.compute_native_head_t(montage)
    
    fig_ = plot_alignment(info, eeg='projected', show_axes=False, 
                          trans=trans, surfaces={"head": 1.0}, coord_frame='mri', 
                          subject=subject_)    
    fig_.plotter.off_screen = True

    set_3d_view(figure=fig_, azimuth=135, elevation=80, distance=0.4)
    fig_.plotter.screenshot(subject_ + "_1.png")

    set_3d_view(figure=fig_, azimuth=45, elevation=80, distance=0.4)
    fig_.plotter.screenshot(subject_ + "_2.png")

    set_3d_view(figure=fig_, azimuth=270, elevation=80, distance=0.4)
    fig_.plotter.screenshot(subject_ + "_3.png")

    
    
    fig_ = plot_alignment(info, eeg=['original', 'projected'], show_axes=True, 
                          trans=trans, surfaces="head", coord_frame='mri', dig=True, 
                          subject=subject_)    
    fig_.plotter.off_screen = True

    set_3d_view(figure=fig_, azimuth=135, elevation=80, distance=0.5, 
                focalpoint=np.linalg.inv(trans["trans"])[:3, 3])
    fig_.plotter.screenshot(subject_ + "_clear_1.png")

    set_3d_view(figure=fig_, azimuth=45, elevation=80, distance=0.5, 
                focalpoint=np.linalg.inv(trans["trans"])[:3, 3])
    fig_.plotter.screenshot(subject_ + "_clear_2.png")

    set_3d_view(figure=fig_, azimuth=270, elevation=80, distance=0.5, 
                focalpoint=np.linalg.inv(trans["trans"])[:3, 3])
    fig_.plotter.screenshot(subject_ + "_clear_3.png")


if __name__ == "__main__":    
    
    Path(path / "images").mkdir(exist_ok=True)
    os.chdir(path / "images")
    for subject in tqdm(subjects):
        print("Processing ", subject, "...")
        montage_fname = fs_subject_path / subject / "HGSN129-montage.fif"
        save_images(subject, str(montage_fname))    
