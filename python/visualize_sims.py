#!/usr/bin/env python

import argparse
import altair as alt
import pandas
import csv
import json
import copy

#import pprint

# Functional helpers; look ye and behold!

def dissoc(d, k):
    d = d.copy()
    try:
        del d[k]
    except:
        pass
    return d

def merge(d, d2):
    d = copy.deepcopy(d)
    d.update(d2)
    return d

def update(d, k, f):
    d = d.copy()
    d[k] = f(d.get(k))
    return d

# Passes through dicts extracted from pandas and converts weird types to simple jsonable things

def clean_data(data):
    if isinstance(data, dict):
        return {k: clean_data(v) for k, v in data.items()}
    elif isinstance(data, list):
        return list(map(clean_data, data))
    elif not pandas.notna(data):
        return None
    else:
        return data



def seqs_chart(posterior_seqs):
    c = alt.Chart(posterior_seqs).mark_bar().encode(
        alt.X('prob:Q', bin=True),
        alt.Y('times_seen:Q', aggregate='sum'),
        alt.Color('truth_dist:N'))
        # for a more semantic blue to red coloring scheme
        #alt.Color('truth_dist:Q', scale=alt.Scale(range=["blue", "red"])))
    c = c.facet(
        row="simulation:N",
        column="seq_type:N")
    return c


def truth_dists_probs_chart(posterior_seqs):
    base = alt.Chart().mark_bar().encode(
        alt.X('prob:Q', bin=True),
        alt.Y('times_seen:Q', aggregate='sum'),
        alt.Color('truth_dist:N'))
        # for a more semantic blue to red coloring scheme
        #alt.Color('truth_dist:Q', scale=alt.Scale(range=["blue", "red"])))
    rule = alt.Chart().mark_rule().encode(
        alt.X('prob:Q', aggregate='mean'),
        alt.Color('truth_dist:N'),
        size=alt.value(3))
    text = alt.Chart().mark_text(color='black').encode(
        alt.X('prob:Q', aggregate='mean'),
        text='truth_dist')
    binned_chart = base + rule + text
    chart = binned_chart.facet(column='seq_type:N') & binned_chart.facet(row='simulation:N', column='seq_type:N')
    chart.data = posterior_seqs
    return chart


def truth_dists_probs_vl(posterior_seqs):
    data = posterior_seqs.to_dict('records')
    main_spec = {'facet': {'column': {'field': 'seq_type', 'type': 'nominal'}},
                 'spec': {'width': 500,
                          'height': 300,
                          'layer': [{'mark': 'bar',
                                     'encoding': {'x': {'field': 'prob',
                                                        'bin': True},
                                                  'y': {'field': 'times_seen',
                                                        'aggregate': "sum"},
                                                  'color': {'field': 'truth_dist'}}},
                                    {'mark': 'rule',
                                     'encoding': {'x': {'field': 'prob',
                                                        'aggregate': 'mean'},
                                                  'color': {'field': 'truth_dist'}}},
                                    {'mark': 'text',
                                     'encoding': {'x': {'field': 'prob',
                                                        'aggregate': 'mean'},
                                                  'text': {'field': 'truth_dist'}}}]}}
    per_tree_spec = update(main_spec, 'facet',
            lambda s: merge(s, {'row': {'field': 'simulation', 'type': 'nominal'}}))
    return {'data': {'values': clean_data(data)},
            'vconcat': [main_spec, per_tree_spec]}



def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('posterior_seqs', type=pandas.read_csv)
    #parser.add_argument('posterior_samples', type=pandas.read_csv)
    parser.add_argument('posterior_samples')
    parser.add_argument('outbase')
    args = parser.parse_args()
    # for now
    args.verbose = True
    return args

def main():
    args = get_args()
    #sc = seqs_chart(args.posterior_seqs)
    #if args.verbose:
        #print("The vega-lite:\n")
        #pprint.pprint(dissoc(sc.to_dict(), 'data'))
    #sc.savechart(args.outbase + '.seqs.html')
    #sc.savechart(args.outbase + '.seqs.json')

    sc = truth_dists_probs_chart(args.posterior_seqs)
    sc.savechart(args.outbase + '.seqs.alt.html')

    #sc = truth_dists_probs_vl(args.posterior_seqs)
    #with open(args.outbase + '.seqs.json', 'w') as fh:
        #json.dump(sc, fh)



if __name__ == '__main__':
    main()


