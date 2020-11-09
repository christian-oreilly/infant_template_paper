#!/home/oreillyc/env_mayada/bin/python

from eegip import Analysis, DatasetMng, SourcesBIDS, Sources

import numpy as np
from pathlib import Path
import os
import mne
import argparse


dataset_name="eegip_london"

subjects = ["ANTS3-0Months3T", "ANTS10-5Months3T", "ANTS9-0Months3T",
            "ANTS12-0Months3T", "ANTS15-0Months3T",
            "ANTS2-0Weeks3T", "ANTS18-0Months3T", "ANTS2-0Years3T",
            "ANTS1-0Months3T", "ANTS2-0Months3T",
            "ANTS4-5Months3T", "ANTS6-0Months3T", "ANTS7-5Months3T"]

config = {"atlas":
           {"vol_labels":
            ['Left-Amygdala',
             'Left-Caudate',
             'Left-Hippocampus',
             'Left-Pallidum',
             'Left-Putamen',
             'Left-Thalamus',
             'Right-Amygdala',
             'Right-Caudate',
             'Right-Hippocampus',
             'Right-Pallidum',
             'Right-Putamen',
             'Right-Thalamus'
             ]
            }
          }

template_check_path = Path("/home/oreillyc/scratch/template_check/")
derivative_root = Path("/home/oreillyc/projects/def-emayada/eegip/london/derivatives/lossless/derivatives/oreillyc")

analysis = Analysis(dataset_name, configs=[config])
montage_name = "HGSN128-montage.fif"


def get_event_types(recording):
    inv_event_id = {val: key for key, val in recording.epochs.event_id.items()}
    return np.array([inv_event_id[id_] for id_ in recording.epochs.events[:, 2]])


def all_files_computed(recording, recording_name, template, event_strs):
    if event_strs is None:
        event_strs = get_event_types(recording)
    for event_str in np.unique(event_strs):

        file_path = template_check_path / "{}_{}_{}_mean.netcdf".format(recording_name, event_str, template)
        if not file_path.exists():
            return False

        file_path = template_check_path / "{}_{}_{}_std.netcdf".format(recording_name, event_str, template)
        if not file_path.exists():
            return False

    return True


def compute_one_smodels(no_subject, recompute=False):
    subject_dir = Path(os.environ["SUBJECTS_DIR"])
    model_type = "average"

    dataset_mng = DatasetMng(config=[config], dataset_name=dataset_name)

    recording_name = [name for name in dataset_mng.dataset.recordings if "m12" not in name][0]

    template = sorted(subjects)[no_subject-1]

    print("Processing template ", template)

    fs_subject = template
    fname_montage = subject_dir / fs_subject / montage_name

    recording = dataset_mng.dataset.recordings[recording_name]

    smodel = recording.subject.add_smodel(name=fs_subject, fs_subject=fs_subject,
                                          type_=model_type, exists="return",
                                          derivative_root=derivative_root)
    smodel.ico = None

    if recompute:
        smodel.recompute()
    else:
        smodel.bem_model
        smodel.bem_solution
        smodel.volume_src
        smodel.surface_src


def compute_smodels(recompute=False):

    dataset_mng = DatasetMng(config=[config], dataset_name=dataset_name)

    for no_template in range(len(subjects)):
        compute_one_smodels(no_template, recompute=recompute)

    print("smodel computed.")
    print("Sources need to be extracted for {} records."
          .format(len([name for name in dataset_mng.dataset.recordings if "m12" not in name])))


def extract_sources_for_a_recording(no_record):
    subject_dir = Path(os.environ["SUBJECTS_DIR"])
    recompute_fmodel = True
    recompute_sources = True
    model_type = "average"
    inv_method = "dSPM"
    snr = 1.0
    lambda2 = None

    dataset_mng = DatasetMng(config=[config], dataset_name=dataset_name)

    recording_name = sorted([name for name in dataset_mng.dataset.recordings if "m12" not in name])[no_record]
    for template in subjects:

        event_strs = ('indirect', 'direct', 'noise')

        print("Processing template ", template)

        fs_subject = template
        derivative_name = "template_validation_" + template

        # The montage is saved with electrodes already placed in the correct positions
        fname_trans = None  # subject_dir / fs_subject / "{}-trans.fif".format(fs_subject)

        recording = dataset_mng.dataset.recordings[recording_name]

        if all_files_computed(recording, recording_name, template, event_strs):
            continue

        recording.montage = mne.channels.read_dig_fif(str(subject_dir / fs_subject / montage_name))
        recording.montage.ch_names = ["E" + str(int(ch_name[3:])) for ch_name in recording.montage.ch_names]
        recording.montage.ch_names[128] = "Cz"

        try:
            recording.get_artifact("preprocessed_raw", recompute=True)
            recording.get_artifact("events", recompute=True)
            recording.get_artifact("epochs", recompute=True)
        except FileNotFoundError:
            continue


        # For cases where one of the event_str has no valid epochs
        if all_files_computed(recording, recording_name, template, event_strs=None):
            continue

        smodel = recording.subject.add_smodel(name=fs_subject, fs_subject=fs_subject,
                                              type_=model_type, exists="return",
                                              derivative_root=derivative_root)

        smodel.ico = None
        smodel.bem_model
        smodel.bem_solution
        smodel.volume_src
        smodel.surface_src        

        print(fs_subject)
        fmodel = recording.add_fmodel(smodel, name=fs_subject, exists="return")
        if fname_trans is not None:
            fmodel.fname_trans = fname_trans
        if recompute_fmodel:
            try:
                fmodel.recompute()
            except FileNotFoundError as e:
                print(e)
                print("Skipping {}. It seems not to be available in lossless qcr version.".format(recording_name))
                continue

        if dataset_mng.dataset_config["is_bids"]:
            sources = SourcesBIDS(fmodel, snr=snr, lambda2=lambda2, method=inv_method,
                                  derivative_name=derivative_name)
        else:
            sources = Sources(fmodel, snr=snr, lambda2=lambda2, method=inv_method)

        if recompute_sources:
            sources.get_data(recompute=True, type_="all")

        # Save full cortex stats
        sources.compute_grouped_sources(grouping="mean")
        sources.compute_grouped_sources(grouping="std")
        # sources.compute_grouped_sources(grouping={"median":np.median})
        # sources.compute_grouped_sources(grouping={"mad":scipy.stats.median_absolute_deviation})
        # For some reason, memory issues were generally occuring while computing the mad. Either,
        # there is a cumulative use of memory from these four aggregate computation or the mad
        # computation is less memory-efficient than the three other statistics.

        # Save label stats
        sources.compute_labels(recompute=recompute_sources, no_epochs=None)

        data = sources.get_data(type_="labels")

        recording.flush()
        fmodel.flush()
        sources.flush()

        data -= data.mean("times")  # De-meaning time series

        event_strs = get_event_types(recording)
        for event_str in np.unique(event_strs):
            epochs_no = np.where(event_strs == event_str)[0]

            file_path = template_check_path / "{}_{}_{}_mean.netcdf".format(recording_name, event_str, template)
            data.sel(epochs=epochs_no).mean("epochs").to_netcdf(file_path, mode='w')

            file_path = template_check_path / "{}_{}_{}_std.netcdf".format(recording_name, event_str, template)
            data.sel(epochs=epochs_no).std("epochs").to_netcdf(file_path, mode='w')

        # Removing full source from the hard drive for space reasons
        sources.delete_sources()


def save_epochs_N():
    dataset_mng = DatasetMng(config=[config], dataset_name=dataset_name)
    recording_names = sorted([name for name in dataset_mng.dataset.recordings if "m12" not in name])

    name_lst = []
    event_lst = []
    N_lst = []
    for recording_name in recording_names:
        recording = dataset_mng.dataset.recordings[recording_name]

        try:
            recording.get_artifact("events", recompute=True)
            recording.get_artifact("epochs", recompute=True)
        except FileNotFoundError:
            continue

        event_strs = get_event_types(recording)

        for event_str in np.unique(event_strs):
            #data.sel(epochs=epochs_no).std("epochs").to_netcdf(file_path, mode='w')
            name_lst.append(recording_name)
            event_lst.append(event_str)
            N_lst.append(np.sum(event_strs == event_str))

    pd.DataFrame({"recording": name_lst, "event": event_lst, "N": N_lst}).to_csv(template_check_path / "N.csv")

def print_report():
    dataset_mng = DatasetMng(config=[config], dataset_name=dataset_name)
    recording_names = sorted([name for name in dataset_mng.dataset.recordings if "m12" not in name])
    for no, recording_name in enumerate(recording_names):
        print(recording_name)
        is_ok = True
        recording = dataset_mng.dataset.recordings[recording_name]
        recording.verbose = False
        try:
            recording.get_artifact("epochs", recompute=False)
        except FileNotFoundError:
            print("{} {}: Skipped. No epochs file".format(no, recording_name))
            continue

        for no_template, template in enumerate(subjects):
            event_strs = ('indirect', 'direct', 'noise')
            if all_files_computed(recording, recording_name, template, event_strs):
                continue

            # For cases where one of the event_str has not valid epochs
            if all_files_computed(recording, recording_name, template, event_strs=None):
                continue

            print("{} {}: {} Missing files for {}".format(no, recording_name, no_template, template))
            is_ok = False

        if is_ok:
            print("{} {} is ok.".format(no, recording_name))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract London06 sources for all infant templates.')
    parser.add_argument('--extract_sources', dest='N', type=int, default=None,
                        help='Recording number from which to extract source. If None, does not extract sources.')
    parser.add_argument('--compute_smodel', dest='compute_smodel', action='store_true',
                        default=False, help='If set, compute the smodel.')
    parser.add_argument('--recompute_smodel', dest='recompute_smodel', action='store_true',
                        default=False, help='Force recomputation of the smodel.')
    parser.add_argument('--recompute_one_smodel', dest='N_model', type=int, default=None,
                        help='Force recomputation of a given smodel.')
    parser.add_argument('--print_report', dest='print_report', action='store_true',
                        default=False, help='Print report.')
    parser.add_argument('--save_epochs_N', dest='save_epochs_N', action='store_true',
                        default=False)

    args = parser.parse_args()
    if args.compute_smodel:
        compute_smodels(False)
    elif args.recompute_smodel:
        compute_smodels(True)
    elif args.print_report:
        print_report()
    elif args.save_epochs_N:
        save_epochs_N()
    else:
        if args.N is not None:
            extract_sources_for_a_recording(args.N)
        if args.N_model is not None:
            compute_one_smodels(args.N_model, True)
