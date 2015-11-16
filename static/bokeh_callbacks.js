function show_data_points(cds, dpds, toggle) {
    var cd = cds.get('data');
    var dpd = dpds.get('data');
    var active = toggle.attributes.active;
    dpd['x'] = [];
    dpd['y'] = [];

    if (active == false) {
        for (var i = 0, i_stop = cd['xs'].length; i < i_stop; i++) {
            dpd['x'] = dpd['x'].concat(cd['xs'][i]);
            dpd['y'] = dpd['y'].concat(cd['ys'][i]);
        }
        toggle.attributes.label = 'Hide data points';
    } else if (active == true) {
        toggle.attributes.label = 'Show data points';
    }
    toggle.trigger('change');
    dpds.trigger('change');
}

function slider(cb_obj, cds, snds, plot_idx_ds, dpds, toggle) {
    var cd = cds.get('data');
    var snd = snds.get('data');
    var dpd = dpds.get('data');
    var current_plot_idx = plot_idx_ds.get('data');
    var plot_idx = cb_obj.get('value');
    var active = toggle.get('active');

    current_plot_idx['i'] = plot_idx;

    for (var key in cd) {
        if (cd.hasOwnProperty(key)) {
            cd[key] = snd['spectral_networks'][plot_idx][key];
        }
    }

    dpd['x'] = [];
    dpd['y'] = [];
    if (toggle.attributes.active == false) {
        for (var i = 0, i_stop = cd['xs'].length; i < i_stop; i++) {
            dpd['x'] = dpd['x'].concat(cd['xs'][i]);
            dpd['y'] = dpd['y'].concat(cd['ys'][i]);
        }
    }

    cds.trigger('change');
    plot_idx_ds.trigger('change');
    dpds.trigger('change');
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
        var Dx = x[a_i] - x[a_i - 1]
        var Dy = y[a_i] - y[a_i - 1]
        cd['arrow_angle'][i] = Math.atan2(Dy, Dx) - (Math.PI / 2);
    }
    cds.trigger('change');
}
