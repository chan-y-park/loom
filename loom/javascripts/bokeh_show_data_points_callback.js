var cd = cds.get('data');
var dpd = dpds.get('data');
dpd['x'] = [];
dpd['y'] = [];

for (var i = 0, i_stop = cd['xs'].length; i < i_stop; i++) {
    dpd['x'] = dpd['x'].concat(cd['xs'][i]);
    dpd['y'] = dpd['y'].concat(cd['ys'][i]);
}
dpds.trigger('change')
