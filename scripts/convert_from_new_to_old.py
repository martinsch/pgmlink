#!/usr/bin/python

import h5py
import numpy as np
import glob


def convert(input_file,
            output_directory,
            output_pattern,
            tg_name,
            label_image_path,
            dtype,
            result_label_path):
    ds = input_file[label_image_path]
    tg = input_file[tg_name]
    for idx, t in enumerate(sorted(tg.keys())):
        print "Processing timestep %s (index %d)" % (t, idx)
        with h5py.File(output_directory + '/' + output_pattern % int(t), 'w') as f:
            f.create_dataset(result_label_path, data=np.squeeze(ds[idx,...]))
            tg_at = f.create_group(tg_name)
            for k in tg[t].keys():
                tg_at.create_dataset(k, data=tg[t][k])

                
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--result-directory', '-r', required=True, help='Path to result directory.')
    parser.add_argument('--result-label-path', '-R', default='segmentation/labels', help='Path to dataset in h5 file. [default=%(default)s]')
    parser.add_argument('--label-image', '-l', default='exported_data', help='Path to label image inside a h5 file. [default=%(default)s]')
    parser.add_argument('--input-file', '-i', required=True, help='Path to input h5 file.')
    parser.add_argument('--dtype', '-t', default='uint32', help='numpy dtype of the result data set. [default=%(default)s]')
    parser.add_argument('--tracking-group', '-g', default='tracking', help='Path to tracking group. [default=%(default)s]')
    parser.add_argument('--output-pattern', '-p', default='%04d.h5', help='Output file pattern. [default=%(default)s]')
    args = parser.parse_args()
    with h5py.File(args.input_file, 'r') as f:
            convert(input_file = f,
                    output_directory = args.result_directory.rstrip('/'),
                    output_pattern = args.output_pattern,
                    tg_name = args.tracking_group,
                    label_image_path = args.label_image,
                    dtype = np.dtype(args.dtype),
                    result_label_path = args.result_label_path)
