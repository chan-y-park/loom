function show_data_points(cds, dpds, hover) {
    var cd = cds.get('data');
    var dpd = dpds.get('data');
    dpd['x'] = [];
    dpd['y'] = [];

    for (var i = 0, i_stop = cd['xs'].length; i < i_stop; i++) {
        dpd['x'] = dpd['x'].concat(cd['xs'][i]);
        dpd['y'] = dpd['y'].concat(cd['ys'][i]);
    }
    hover.attributes.tooltips = null;
    dpds.trigger('change');
}

function hide_data_points(cds, dpds, hover) {
    var cd = cds.get('data');
    var dpd = dpds.get('data');
    dpd['x'] = [];
    dpd['y'] = [];

    hover.attributes.tooltips = [['name', '@label'], ['root', '@root']];
    dpds.trigger('change');
}

function sn_slider(
    cb_obj, cds, snds, sn_idx_ds, dpds, pds, hover,
    plot_options_ds, tree_idx_ds
) {
    var cd = cds.get('data');
    var snd = snds.get('data');
    var dpd = dpds.get('data');
    var pd = pds.get('data');
    var current_sn_idx = sn_idx_ds.get('data');
    var sn_idx = cb_obj.get('value');
    var plot_options = plot_options_ds.get('data');
    var notebook = plot_options['notebook'];
    var show_trees = plot_options['show_trees'];
    var tree_idx = tree_idx_ds.get('data');
    
    current_sn_idx['i'] = sn_idx;
    tree_idx['j'] = 0;
    if (show_trees == 'true') {
        document.getElementById("current_tree_idx").innerHTML = 'All';
    }

    for (var key in cd) {
        if (cd.hasOwnProperty(key)) {
            if (show_trees == 'false') {
                cd[key] = snd['spectral_networks'][sn_idx][key];
            } else {
                cd[key] = snd['spectral_networks'][sn_idx][0][key];
            }
        }
    }
    cds.trigger('change');
    sn_idx_ds.trigger('change');
    tree_idx_ds.trigger('change');
    hide_data_points(cds, dpds, hover);
    if (notebook == 'false') {
        document.getElementById("phase").innerHTML = pd['phase'][sn_idx];
    }
}

function change_soliton_tree(
    cds, snds, sn_idx_ds, tree_idx_ds, plot_options_ds, change
) {
    var cd = cds.get('data');
    var snd = snds.get('data');
    var sn_idx = sn_idx_ds.get('data');
    var tree_idx = tree_idx_ds.get('data');
    var plot_options = plot_options_ds.get('data');
    var notebook = plot_options['notebook'];
    var show_trees = plot_options['show_trees'];

    var sn_i = sn_idx['i'];
    var tree_j = Number(tree_idx['j']) + change;
    var max_tree_j = snd['spectral_networks'][sn_i].length;
    if (tree_j > max_tree_j) {
        tree_j = max_tree_j;
    } else if (tree_j < 0) {
        tree_j = 0;
    }
    
    for (var key in cd) {
        if (cd.hasOwnProperty(key)) {
            cd[key] = snd['spectral_networks'][sn_i][tree_j][key];
        }
    }

    if (notebook == 'false') {
        var tree_idx_label = document.getElementById("current_tree_idx");
        if (tree_j == 0) {
            tree_idx_label.innerHTML = 'All';
        } else {
            tree_idx_label.innerHTML = tree_j - 1;
        }
    }
    
    cds.trigger('change');
    tree_idx['j'] = String(tree_j);
    tree_idx_ds.trigger('change');
}

function show_prev_soliton_tree(
    cds, snds, sn_idx_ds, tree_idx_ds, plot_options_ds
) {
    change_soliton_tree(
        cds, snds, sn_idx_ds, tree_idx_ds, plot_options_ds, -1
    );
}

function show_next_soliton_tree(
    cds, snds, sn_idx_ds, tree_idx_ds, plot_options_ds
) {
    change_soliton_tree(
        cds, snds, sn_idx_ds, tree_idx_ds, plot_options_ds, 1
    );   
}

function redraw_arrows(cds, x_range, y_range) {
    var cd = cds.get('data');
    var x_s = x_range.get('start');
    var x_e = x_range.get('end');
    var y_s = y_range.get('start');
    var y_e = y_range.get('end');

    for (var i = 0, i_stop = cd['arrow_x'].length; i < i_stop; i++) {
        // Domain of the segment.
        var range = cd['ranges'][i];
        var x_min = range[0];
        var x_max = range[1];
        var y_min = range[2];
        var y_max = range[3];

        if ((x_max < x_s) || (x_min > x_e) || (y_max < y_s) || (y_min > y_e)) {
            // The segment is out of screen.
            continue;
        }

        // Now find the new location for the arrow.
        var a_x;
        var a_y;
        var a_i;
        var x = cd['xs'][i];
        var y = cd['ys'][i];
        var x_length = x.length;
        var denom = 2;
        var found = false;

        while ((Math.floor(x_length / denom) > 0) && !found) {
            for (var j = 1; j < denom; j += 2) {
                a_i = Math.floor(x_length * j / denom);
                a_x = x[a_i];
                a_y = y[a_i];
                if ((a_x >= x_s) && (a_x <= x_e) && 
                    (a_y >= y_s) && (a_y <= y_e)) {
                    found = true;
                    break;
                }
            }
            denom *= 2;
        }
        cd['arrow_x'][i] = a_x;
        cd['arrow_y'][i] = a_y;
        var Dx = x[a_i] - x[a_i - 1];
        var Dy = y[a_i] - y[a_i - 1];
        cd['arrow_angle'][i] = Math.atan2(Dy, Dx) - (Math.PI / 2);
    }
    cds.trigger('change');
}

function update_plot_range(x_range, y_range) {
    var x_s = x_range.get('start');
    var x_e = x_range.get('end');
    var y_s = y_range.get('start');
    var y_e = y_range.get('end');

    var plot_range_input = document.getElementById("plot_range_input");
    plot_range_input.value =
        "[[" + x_s + "," + x_e + "],[" + y_s + "," + y_e + "]]";
}
