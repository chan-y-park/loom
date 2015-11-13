var cd = cds.get('data');
var snd = snds.get('data');
var current_plot_idx = plot_idx_ds.get('data');
var plot_idx = cb_obj.get('value');

current_plot_idx['i'] = plot_idx;

for (var key in cd) {
    if (cd.hasOwnProperty(key)) {
        cd[key] = snd['spectral_networks'][plot_idx][key];
    }
}

cds.trigger('change');
plot_idx_ds.trigger('change');
