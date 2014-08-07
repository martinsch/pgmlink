import pgmlink
import numpy as np
import h5py
import os
import vigra

def write_dict_value(dic, key, value):
    if len(value) == 0:
        return
    else:
        dic[key] = value
    return dic


def get_events(eventsVector):
    events = {}
    for t in range(len(eventsVector)):
        events[str(t)] = get_events_at(eventsVector, t)
    return events

def get_events_at(eventsVector, t):  
    dis = []
    app = []
    div = []
    mov = []
    merger = []
                
    for event in eventsVector[t]:
        if event.type == pgmlink.EventType.Appearance:
            app.append((event.traxel_ids[0], event.energy))
        if event.type == pgmlink.EventType.Disappearance:
            dis.append((event.traxel_ids[0], event.energy))
        if event.type == pgmlink.EventType.Division:
            div.append((event.traxel_ids[0], event.traxel_ids[1], event.traxel_ids[2], event.energy))
        if event.type == pgmlink.EventType.Move:
            mov.append((event.traxel_ids[0], event.traxel_ids[1], event.energy))
        if hasattr(pgmlink.EventType, "Merger") and event.type == pgmlink.EventType.Merger:                    
            merger.append((event.traxel_ids[0], event.traxel_ids[1], event.energy))

    # convert to ndarray for better indexing
    events_at = {}
    write_dict_value(events_at, "dis", np.asarray(dis))
    write_dict_value(events_at, "app", np.asarray(app))
    write_dict_value(events_at, "div", np.asarray(div))
    write_dict_value(events_at, "mov", np.asarray(mov))
    write_dict_value(events_at, "merger", np.asarray(merger))

    return events_at


def dump_tif(image, out_dir, iteration_axis = -1):
    n_images = image.shape[iteration_axis]
    print "Dumping image with shape", image.shape, "to %d images in directory %s" % (n_images, out_dir) 
    n_digits = len(str(n_images))
    base_string = out_dir.rstrip('/') + '/' + 'exported_%' + ('0%d' % n_digits) + 'd.tif'
    slicing = [slice(None)] * len(image.shape)
    for i in range(n_images):
        print "Writing " + base_string % i
        slicing[iteration_axis] = i
        vigra.impex.writeImage(image[tuple(slicing)]*(image[tuple(slicing)] > 0), base_string % i, dtype='UINT16')


def project_active_file(fn, LEVEL_AXIS=-1, TIME_AXIS=-2, timestep_begin=0, timestep_end=-1, level_path='exported_data'):
    ds_name = 'projected_active_regions'
    print "projecting active regions in file %s/%s" % (fn, ds_name)
    with h5py.File(fn, 'r+') as f:
        shape = list(f[level_path].shape)
        if timestep_begin < 0:
            timestep_begin = shape[TIME_AXIS] + timestep_begin
        if timestep_end < 0:
            timestep_end = shape[TIME_AXIS] + timestep_end + 1
        n_timesteps = timestep_end - timestep_begin
        del shape[LEVEL_AXIS]
        if ds_name in f.keys():
            del f[ds_name]
        active_regions = f.create_dataset(ds_name, shape = shape, dtype = np.int32)
        if len(shape) == 4:
            axistags = 'txyz'
        else:
            axistags = 'txy'
        active_regions.attrs['axistags'] = vigra.defaultAxistags(axistags).toJSON()
        project_active(f[level_path], active_regions, f['tracking'], LEVEL_AXIS, TIME_AXIS, timestep_begin, timestep_end)


def project_active(level_image, active_regions, tracking_group, LEVEL_AXIS, TIME_AXIS, timestep_begin, timestep_end):
    shape = list(level_image.shape)
    if LEVEL_AXIS < 0:
        LEVEL_AXIS = len(shape) + LEVEL_AXIS
    if TIME_AXIS < 0:
        TIME_AXIS = len(shape) + TIME_AXIS
    if timestep_begin < 0:
        timestep_begin = shape[TIME_AXIS] + timestep_begin
    if timestep_end < 0:
        timestep_end = shape[TIME_AXIS] + timestep_end + 1
    transfer_labels = get_label_dict(tracking_group)
    for idx, ts in enumerate(xrange(timestep_begin, timestep_end)):
        print "projecting timestep %d (index %d):" % (ts, idx)
        slices = [slice(None)] * len(level_image.shape)
        slices[TIME_AXIS] = ts
        image_at = level_image[tuple(slices)]
        unique_labels = np.arange(0, np.amax(level_image) + 1, dtype = np.int32)
        unique_labels[0] = 0
        for i, val in enumerate(unique_labels):
            if i == 0:
                continue
            if i not in transfer_labels[ts].keys():
                unique_labels[i] = 0

        active_regions_slices = [slice(None)] * len(active_regions.shape)
        active_regions_slices[TIME_AXIS] = ts
        level_axis = LEVEL_AXIS
        if TIME_AXIS < LEVEL_AXIS:
            level_axis -= 1
        active_regions[tuple(active_regions_slices)] = unique_labels[image_at].max(axis=level_axis)
            
        
    

def get_number_of_timesteps(shape, TIME_AXIS, ts_first, ts_last):
    if ts_first < 0:
        ts_first = shape[TIME_AXIS] + ts_first
    if ts_last < 0:
        ts_last = shape[TIME_AXIS] + ts_last + 1
    assert(ts_last >= ts_first)
    return ts_last - ts_first
    
def relabel_file(fn, LEVEL_AXIS=-1, TIME_AXIS=-2, n_timesteps=None):
    print "relabeling data in file %s" % fn
    with h5py.File(fn, 'r+') as f:
        shape = list(f['exported_data'].shape)
        if n_timesteps is not None:
            shape[TIME_AXIS] = n_timesteps
        if "lineage_with_channels" in f.keys():
            del f['lineage_with_channels']
        lineage_with_channels = f.create_dataset('lineage_with_channels', shape = f['exported_data'].shape, dtype=np.int64)
        del shape[LEVEL_AXIS]
        if "lineage" in f.keys():
            del f['lineage']
        lineage = f.create_dataset("lineage", shape=shape, dtype=np.int64)
        # relabel_deprecated(f['exported_data'], lineage, lineage_with_channels, f['tracking'], LEVEL_AXIS, TIME_AXIS)
        relabel(f['exported_data'], lineage, f['tracking'], LEVEL_AXIS, TIME_AXIS, n_timesteps)
        # temporary . as output path
        # dump_tif(lineage, '.', )


def get_label_dict(tracking_group):
    transfer_labels = []
    # transfer_labels.append({})
    keys_sorted = sorted([int(x) for x in tracking_group.keys()])
    for idx, ts in enumerate(keys_sorted[:]):
        transfer_labels.append({})
        tg_at = tracking_group['%04d' % ts]
        if "Appearances" in tg_at.keys():
            for app in tg_at["Appearances"]:
                transfer_labels[-1][app[0]] = np.random.randint(1, 255)
                
        if "Moves" in tg_at.keys():
            for mov in tg_at["Moves"]:
                if not transfer_labels[-2].has_key(mov[0]):
                    transfer_labels[-2][mov[0]] = np.random.randint(1, 255)
                transfer_labels[-1][mov[1]] = transfer_labels[-2][mov[0]]

        if "Splits" in tg_at.keys():
            for div in tg_at["Splits"]:
                if not transfer_labels[-2].has_key(div[0]):
                    transfer_labels[-2][div[0]] = np.random.randint(1, 255)
                transfer_labels[-1][div[1]] = transfer_labels[-2][div[0]]
                transfer_labels[-1][div[2]] = transfer_labels[-2][div[0]]

    return transfer_labels

    
def relabel(level_image, label_image, tracking_group, LEVEL_AXIS=-1, TIME_AXIS=-2, n_timesteps=None):
    transfer_labels = get_label_dict(tracking_group)
    if n_timesteps is None:
        n_timesteps = level_image.shape[TIME_AXIS]
    for ts in xrange(n_timesteps):
        print "relabeling timestep %d:" % ts
        slices = [slice(None)] * len(level_image.shape)
        slices[TIME_AXIS] = ts
        image_at = level_image[tuple(slices)]
        unique_labels = np.arange(0, np.amax(level_image) + 1, dtype = np.int64)
        unique_labels[0] = -2
        unique_labels[1:] = -1
        unique_labels[transfer_labels[ts].keys()] = transfer_labels[ts].values()
        label_image_slices = [slice(None)] * len(label_image.shape)
        label_image_slices[TIME_AXIS] = ts
        label_image[tuple(label_image_slices)] = unique_labels[image_at].max(axis=LEVEL_AXIS)
        
        

def relabel_deprecated(level_image, label_image, label_image_with_levels, tracking_group, LEVEL_AXIS=-1, TIME_AXIS=-2):
    max_label = 1
    transfer_labels = []
    keys_sorted = sorted([int(x) for x in tracking_group.keys()])
    assert len(keys_sorted) <= level_image.shape[-2]
    for ts, k in enumerate(keys_sorted):
        print "relabeling timestep %d (index %d):" % (k, ts)
        k_str = str(k)
        slices = [slice(None)] * len(level_image.shape)
        slices[TIME_AXIS] = k
        image_at = level_image[tuple(slices)]
        image_at = np.require(image_at, np.int64)
        print "image dimension -", image_at.shape
        unique_labels = np.arange(0, np.amax(image_at)+1, dtype=image_at.dtype)
        unique_labels[1:] = -1
        unique_labels[0] = -2

        tg_at = tracking_group[k_str]
        
        app = []
        dis = []
        mov = []
        div = []
        if "Appearances" in tg_at.keys():
            app = tg_at["Appearances"][...]
        if "Disappearances" in tg_at.keys():
            dis = tg_at["Disappearances"][...]
        if "Moves" in tg_at.keys():
            mov = tg_at["Moves"][...]
        if "Splits" in tg_at.keys():
            div = tg_at["Splits"][...]

        transfer_labels.append({})
        curr = transfer_labels[-1]
        assert np.all(np.array(transfer_labels[-1].values()[1:]) == -1), transfer_labels[-1].values()[1:]

        if k == keys_sorted[0]:
            print "Initializing positive detections in the earliest timestep"
            # negative detections are set to -1 by default anyway
            for e in mov:
                transfer_labels[-1][e[0]] = e[0]
            for e in dis:
                transfer_labels[-1][e[0]] = e[0]
            for e in div:
                transfer_labels[-1][e[0]] = e[0]
            max_label = max(max_label, unique_labels.shape[0])+1
            
        else:
            tg_prev = tracking_group[str(keys_sorted[ts-1])]
            prev = transfer_labels[-2]
            count = 0
            for e in app:
                # print "Initializing appearance at time %d with id %d to lineage id %d" % (ts, e[0], max_label)
                transfer_labels[-1][e[0]] = max_label
                max_label += 1
                count += 1
            if "Moves" in tg_prev.keys():
                for e in tg_prev["Moves"]:
                    # print "Move from %d@%d to %d@%d on lineage %d" % (e[0], ts-1, e[1], ts, transfer_labels[-2][e[0]])
                    assert transfer_labels[-2][e[0]] > 0
                    transfer_labels[-1][e[1]] = transfer_labels[-2][e[0]]
                    assert transfer_labels[-1][e[1]] > 0
                    count += 1
                    # print "mov:", e, transfer_labels[-2][e[0]], transfer_labels[-1][e[1]]
            if "Splits" in tg_prev.keys():
                for e in tg_prev["Splits"]:
                    # print "Division from %d@%d to (%d,%d)@%d on lienage %d" % (e[0], ts-1, e[1], e[2], ts, transfer_labels[-2][e[0]])
                    assert transfer_labels[-2][e[0]] > 0
                    transfer_labels[-1][e[1]] = transfer_labels[-2][e[0]]
                    transfer_labels[-1][e[2]] = transfer_labels[-2][e[0]]
                    assert transfer_labels[-1][e[1]] > 0
                    assert transfer_labels[-1][e[2] ]> 0
                    # print "div:", e, transfer_labels[-2][e[0]], transfer_labels[-1][e[1]], transfer_labels[-1][e[2]]
                    count += 2

        unique_labels[transfer_labels[-1].keys()] = transfer_labels[-1].values()
        relabeled_image = unique_labels[image_at]
        label_image_with_levels[tuple(slices)] = relabeled_image
        label_image[..., k] = relabeled_image.max(axis=-1)

        assert len(transfer_labels) == ts+1


def colorize_file(fn, TIME_AXIS=-1):
    with h5py.File(fn, 'r+') as f:
        if "colorized" in f.keys():
            del f['colorized']
        colorized = f.create_dataset('colorized', data=np.zeros(f['lineage'].shape+(3,)), dtype=np.uint8)
        colorize(f['lineage'], colorized, TIME_AXIS)

        
def colorize_deprecated(lineage, colorized, TIME_AXIS=-1):
    timesteps = lineage.shape[TIME_AXIS]
    colors = {}
    colors[0] = np.array([0, 0, 0], dtype=np.uint8)
    colors[-1] = np.array([50, 50, 50], dtype=np.uint8)
    for t in xrange(timesteps):
        print "Coloring timestep %d" % t
        slices = [slice(None)] * len(lineage.shape)
        slices[TIME_AXIS] = t
        slices_with_channels = slices + [slice(None)]
        image_at = lineage[tuple(slices)]
        colored_at = np.zeros(image_at.shape + (3,))
        labels = np.arange(0, np.amax(image_at)+1, dtype=image_at.dtype)
        for l in labels:
            if l == 0:
                continue
            if l not in colors.keys():
                colors[l] = np.random.random_integers(0, 255, (3,)).astype(np.uint8)
        

        color_arr = np.zeros((len(colors.keys()), 3), dtype=np.uint8)
        color_arr[colors.keys(),:] = colors.values()
        colored_at = color_arr[image_at]
        colorized[tuple(slices_with_channels)] = colored_at
            
            
def non_zero_lineage_file(fn, iteration_axis=-1):
    with h5py.File(fn, 'r+') as f:
        if 'lineage' not in f:
            print 'WARNING: will not write lineage_non_zero since there is no lineage'
            return
        if 'lineage_non_zero' in f:
            del f['lineage_non_zero']
        lineage = f['lineage']
        lineage_non_zero = f.create_dataset('lineage_non_zero', dtype = np.uint32, shape = lineage.shape)
        non_zero_lineage(lineage, lineage_non_zero, iteration_axis)


def non_zero_lineage(lineage, lineage_non_zero, iteration_axis=-1):
    slicing = [slice(None)] * len(lineage.shape)
    for i in xrange(lineage.shape[iteration_axis]):
        slicing[iteration_axis] = i
        lineage_at = lineage[tuple(slicing)]
        lineage_non_zero[tuple(slicing)] = lineage_at*(lineage_at >= 0)
                
            
            

def write_events(events, fn):
    
    print "-- Writing results to " + os.path.basename(fn) + '/tracking'
    # write only if file exists
    with h5py.File(fn, 'r+') as f_curr:
        # delete old tracking
        if "tracking" in f_curr.keys():
            del f_curr["tracking"]
        g = f_curr.create_group("tracking")
            
        for ts, events_at in enumerate(events):
            print "  -- Writing results for timestep % 4d" % ts
            tg = g.create_group("%04d" % ts)
            dis = []
            app = []
            div = []
            mov = []
            mer = []
            mul = []
            for event in events_at:
                if event.type == pgmlink.EventType.Appearance:
                    app.append((event.traxel_ids[0], event.energy))
                if event.type == pgmlink.EventType.Disappearance:
                    dis.append((event.traxel_ids[0], event.energy))
                if event.type == pgmlink.EventType.Division:
                    div.append((event.traxel_ids[0], event.traxel_ids[1], event.traxel_ids[2], event.energy))
                if event.type == pgmlink.EventType.Move:
                    mov.append((event.traxel_ids[0], event.traxel_ids[1], event.energy))
                if event.type == pgmlink.EventType.Merger:
                    mer.append((event.traxel_ids[0], event.traxel_ids[1], event.energy))
                if event.type == pgmlink.EventType.MultiFrameMove:
                    mul.append(tuple(event.traxel_ids) + (event.energy,))

            # convert to ndarray for better indexing
            dis = np.asarray(sorted(dis, key=lambda el: el[0]))
            app = np.asarray(sorted(app, key=lambda el: el[0]))
            div = np.asarray(sorted(div, key=lambda el: el[0]))
            mov = np.asarray(sorted(mov, key=lambda el: el[0]))
            mer = np.asarray(sorted(mer, key=lambda el: el[0]))
            mul = np.asarray(sorted(mul, key=lambda el: el[0]))
    
        

        
        
            # write associations
            if len(app):
                ds = tg.create_dataset("Appearances", data=app[:, :-1], dtype=np.int32)
                ds.attrs["Format"] = "cell label appeared in current file"
            
                ds = tg.create_dataset("Appearances-Energy", data=app[:, -1], dtype=np.double)
                ds.attrs["Format"] = "lower energy -> higher confidence"

            if len(dis):
                ds = tg.create_dataset("Disappearances", data=dis[:, :-1], dtype=np.int32)
                ds.attrs["Format"] = "cell label disappeared in current file"
                
                ds = tg.create_dataset("Disappearances-Energy", data=dis[:, -1], dtype=np.double)
                ds.attrs["Format"] = "lower energy -> higher confidence"


            if len(mov):
                ds = tg.create_dataset("Moves", data=mov[:, :-1], dtype=np.int32)
                ds.attrs["Format"] = "from (previous file), to (current file)"
            
                ds = tg.create_dataset("Moves-Energy", data=mov[:, -1], dtype=np.double)
                ds.attrs["Format"] = "lower energy -> higher confidence"
        
                
            if len(div):
                ds = tg.create_dataset("Splits", data=div[:, :-1], dtype=np.int32)
                ds.attrs["Format"] = "ancestor (previous file), descendant (current file), descendant (current file)"
            
                ds = tg.create_dataset("Splits-Energy", data=div[:, -1], dtype=np.double)
                ds.attrs["Format"] = "lower energy -> higher confidence"

            if len(mer):
                ds = tg.create_dataset("Mergers", data=mer[:, :-1], dtype=np.int32)
                ds.attrs["Format"] = "descendant (current file), number of objects"
		
                ds = tg.create_dataset("Mergers-Energy", data=mer[:, -1], dtype=np.double)
                ds.attrs["Format"] = "lower energy -> higher confidence"

            if len(mul):
                ds = tg.create_dataset("MultiFrameMoves", data=mul[:, :-1], dtype=np.int32)
                ds.attrs["Format"] = "from (given by timestep), to (current file), timestep"
                    
                ds = tg.create_dataset("MultiFrameMoves-Energy", data=mul[:, -1], dtype=np.double)
                ds.attrs["Format"] = "lower energy -> higher confidence"
                

    print "-> results successfully written"


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--filename', '-f', required=True, help='Path to some h5 file.')
    parser.add_argument('--dataset-path', '-d', default='exported_data', help='Internal path to dataset [default=].')
    parser.add_argument('--level-axis', '-l', type=int, default=-1, help='Level axis [default=].')
    parser.add_argument('--time-axis', '-T', type=int, default=0, help='Time axis [default=].')
    parser.add_argument('--timestep-begin', '-b', type=int, default=0, help='First timestep [default=].')
    parser.add_argument('--timestep-end', '-e', type=int, default=-1, help='Last timestep [default=].')
    options = parser.parse_args()

    project_active_file(fn=options.filename,
                        LEVEL_AXIS=options.level_axis,
                        TIME_AXIS=options.time_axis,
                        timestep_begin=options.timestep_begin,
                        timestep_end=options.timestep_end,
                        level_path=options.dataset_path)
