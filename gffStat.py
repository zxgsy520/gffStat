#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import sys
import copy

#print(sys.path)
import logging
import argparse
from GffReader import open_gff
from plotDistribution import plotDistribution


LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("fanjunpeng, Xingguo Zhang",)
__email__ = "jpfan@whu.edu.cn, invicoun@foxmail.com"
__all__ = []


class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def gff2dict(fn):
    """

    :param fn:
    :return:
    """

    # link gene to mRNA, mRNA to CDS, exon by parent
    r = {}
    gene_level = {}
    mRNA_level = {}

    for gff in open_gff(fn):
        _type = gff.type
        _attrs = gff.attributes
        _id = None
        _parents = None
        locus_tag = None
        

        if "ID" in _attrs:
            _id = _attrs["ID"]
        if "Parent" in _attrs:
            _parents = _attrs["Parent"]

        if _type == "gene":
            if not _id:
                print("gene %r has no ID, pass" % gff.to_string())
                continue
            if _id in r:  # gene named after mRNA....
                r[_id]["self"] = gff
                continue

            gene_level[_id] = {"self": gff, "mRNA": []}  # construct gene

        if _type in "mRNA":
            if not _id:
                print("mRNA %r has no ID, pass" % gff.to_string())
                continue

            if not _parents:  # mRNA has no parent, reconstruct gene
                # print("mRNA %r has no parent, reconstruct" % _id)
                gene = copy.deepcopy(gff)
                gid = _id+"-gene"
                gene.type = "gene"
                gene.attributes = {"ID": gid}
                gene_level[gid] = {"self": gene, "mRNA": [_id]}

                _attrs["Parent"] = gid
                mRNA_level[_id] = {"self": gff, "CDS": [], "exon": []}
                continue

            if "," in _parents:
                print("mRNA %r belongs to muti gene! pass it" % _id)
                continue

            if _parents not in gene_level:  # name gene by mRNA parent
                gene_level[_parents] = {"self": None, "mRNA": [_id]}
            else:
                gene_level[_parents]["mRNA"].append(_id)

            if _id in mRNA_level:
                mRNA_level[_id]["self"] = gff
            else:
                mRNA_level[_id] = {"self": gff, "CDS": [], "exon": []}

        if _type in "CDS|exon":
            if not _parents:
                print("%s '%s' has no parent, pass" % (_type, gff.to_string()))
                continue

            for parent in _parents.split(","):
                if parent not in mRNA_level:  # name mRNA by CDS parent
                    mRNA_level[parent] = {"self": None, "CDS": [], "exon": []}
                    mRNA_level[parent][_type].append(gff)
                    continue

                mRNA_level[parent][_type].append(gff)

    return gene_level, mRNA_level


def default_filter(gene_level, mRNA_level):

    filtered_gene_level = {}
    filtered_mRNA_level = mRNA_level

    for gid in gene_level:
        mids = gene_level[gid]["mRNA"]

        if not mids:
            print("%s has no mRNA, pass" % gid)
            continue

        mrna_status = 0
        for mid in mids:

            if mid not in mRNA_level:
                print(("%s has no mRNA, pass" % gid))
                gene_level[gid]["mRNA"].remove(mid)
                continue

            if not mRNA_level[mid]["CDS"]:
                print(("%s has no CDS, pass" % gid))
                gene_level[gid]["mRNA"].remove(mid)
                continue
            if not mRNA_level[mid]["CDS"]:
                print(("%s has no record, pass" % gid))
                gene_level[gid]["mRNA"].remove(mid)
                continue

            mrna_status = 1

        if mrna_status != 0:
            filtered_gene_level[gid] = gene_level[gid]

    return filtered_gene_level, filtered_mRNA_level


def longest_mRNA_filter(gene_level, mRNA_level):
    """

    :param gene_level:
    :param mRNA_level:
    :return:
    """
    print("Filter gff3 by longest mRNA")

    filtered_gene_level = {}
    filtered_mRNA_level = mRNA_level

    for gid in gene_level:
        mids = gene_level[gid]["mRNA"]

        filtered_gene_level[gid] = gene_level[gid]
        longest_mRNA = None

        for mid in mids:

            if not longest_mRNA:
                longest_mRNA = mid
            elif mRNA_level[mid]["self"].length > mRNA_level[longest_mRNA]["self"].length:
                longest_mRNA = mid
            else:
                pass

        filtered_gene_level[gid]["mRNA"] = [longest_mRNA]
        try:
            filtered_gene_level[gid]["self"].start = mRNA_level[longest_mRNA]["self"].start
            filtered_gene_level[gid]["self"].end = mRNA_level[longest_mRNA]["self"].end
        except:
            print("%s ?" % gid)


    return filtered_gene_level, filtered_mRNA_level


def no_UTR_filter(gene_level, mRNA_level):
    """

    :param gene_level:
    :param mRNA_level:
    :return:
    """
    print("Filter gff3 by no UTR")

    filtered_gene_level = {}
    filtered_mRNA_level = {}

    for mid, mv in mRNA_level.items():

        cds = sorted(mv["CDS"], key=lambda i: i.start)
        try:
            mv["self"].start = cds[0].start
            mv["self"].end = cds[-1].end
        except:
            print("%r is not a mRNA, maybe a ncRNA." % mid)
            continue
        mv["CDS"] = cds
        mv["exon"] = []
        filtered_mRNA_level[mid] = mv

    for gid, gv in gene_level.items():
        mids = gene_level[gid]["mRNA"]
        if not mids:
            print("%s has no mRNA, pass" % gid)
            continue

        start = end = 0
        for mid in mids:

            if start == 0 or end == 0:
                start = filtered_mRNA_level[mid]["self"].start
                end = filtered_mRNA_level[mid]["self"].end
                continue

            if filtered_mRNA_level[mid]["self"].start < start:
                start = filtered_mRNA_level[mid]["self"].start
            if filtered_mRNA_level[mid]["self"].end > end:
                end = filtered_mRNA_level[mid]["self"].end

        if hasattr(gv["self"], "start"):
            gv["self"].start = start
            gv["self"].end = end
            filtered_gene_level[gid] = gv
        else:
            logging.warning("%s format error" % gid)
            

    return filtered_gene_level, filtered_mRNA_level


def gff_format(gene_level, mRNA_level):
    r = {}

    # reconstruct exon and intron if necessary
    new_mRNA_level = {}
    for mid, mv in mRNA_level.items():
        if not mv["self"]: # mRNA has no named, pass
            continue

        if not mv["CDS"]:  # mRNA has no CDS, pass
            continue

        new_mRNA_level[mid] = mv
        cds = sorted(mv["CDS"], key=lambda i: i.start)
        new_mRNA_level[mid]["CDS"] = cds

        n = len(cds)
        if not mv["exon"]:  # construct exon
            # print("mRNA %s has no exon, reconstruct from CDS and mRNA" % mid)

            for i in range(n):
                exon = copy.deepcopy(cds[i])
                exon.type = "exon"
                exon.phase = "."
                if i == 0:
                    exon.start = mv["self"].start
                if i == n - 1:
                    exon.end = mv["self"].end
                new_mRNA_level[mid]["exon"].append(exon)
        else:
            new_mRNA_level[mid]["exon"] = sorted(mv["exon"], key=lambda i: i.start)

        # construct intron
        new_mRNA_level[mid]["intron"] = []
        for i in range(n):
            if i == 0:
                continue
            start = cds[i - 1].end + 1
            end = cds[i].start - 1
            if start >= end:
                continue

            intron = copy.deepcopy(cds[i])
            intron.type = "intron"
            intron.phase = "."
            intron.start = start
            intron.end = end
            if intron.length <= 0:
                print("intron %s <= 0" % intron.to_string())
                continue
            new_mRNA_level[mid]["intron"].append(intron)

    # add mRNA to gene
    for gid, gv in gene_level.items():
        if not gv["mRNA"]:
            print("%s has no mRNA, pass" % gid)
            continue

        if not gv["self"]:
            print("%s has no gff, pass" % gid)
            continue

        r[gid] = {"self": gv["self"], "mRNA": {}}
        for mid in gv["mRNA"]:
            r[gid]["mRNA"][mid] = new_mRNA_level[mid]
   # print(r)
    return r


def gff2str(gffdict, fn):

    out = open(fn, "w")
    for gid, gv in sorted(gffdict.items(), key=lambda d:(d[1]["self"].seqid, d[1]["self"].start)):
        out.write(gv["self"].to_string()+"\n")
        for mid, mv in gv["mRNA"].items():
            out.write(mv["self"].to_string()+"\n")
            for exon in mv["exon"]:
                out.write(exon.to_string()+"\n")
            for cds in mv["CDS"]:
                out.write(cds.to_string()+"\n")
    out.close()

    return fn


def stat(gffdict):
    r = {"gene_num": len(gffdict), "gene_length": [],
         "CDS_num": [], "CDS_length": [],
         "intron_num": [], "intron_length": [],
         "exon_num": [], "exon_length": [], }

    for gid, gv in sorted(gffdict.items(), key=lambda d:(d[1]["self"].seqid, d[1]["self"].start)):
        r["gene_length"].append(gv["self"].end-gv["self"].start+1)

        for mid, mv in gv["mRNA"].items():
            r["CDS_num"].append(len(mv["CDS"]))
            CDS_length = sum([i.end-i.start+1 for i in mv["CDS"]])
            r["CDS_length"].append(CDS_length)
            if mv["intron"]:
                r["intron_num"].append(len(mv["intron"]))
                intron_length = [i.end-i.start+1 for i in mv["intron"]]
                r["intron_length"] += intron_length
            else:
                r["intron_num"].append(0)

            r["exon_num"].append(len(mv["exon"]))
            exon_length = [i.end-i.start+1 for i in mv["exon"]]
            r["exon_length"] += exon_length

    return r


def gff_stat(fn, outdir):

    print("Loading gff3 information from %s" % fn)
    gene, mRNA = gff2dict(fn)

    filtered_gene_level, filtered_mRNA_level = default_filter(gene, mRNA)

    filtered_gene_level, filtered_mRNA_level = longest_mRNA_filter(filtered_gene_level, filtered_mRNA_level)
    filtered_gene_level, filtered_mRNA_level = no_UTR_filter(filtered_gene_level, filtered_mRNA_level)
    print("Format gff3 loaded")
    filtered_gff = gff_format(filtered_gene_level, filtered_mRNA_level)
    print("Statistics on gff3")
    stats = stat(filtered_gff)
    print("Output filted gff3")
    gff2str(filtered_gff, "%s/%s.longest.gff3" % (outdir, os.path.basename(fn)))
    return stats


def stat2str(stats):
    if stats["gene_num"]:
        mean_gene_length = sum(stats["gene_length"])/1.0/len(stats["gene_length"])
    else:
        mean_gene_length = 0

    if stats["CDS_num"]:
        mean_cds_num = sum(stats["CDS_num"])/1.0/len(stats["CDS_num"])
        mean_cds_length = sum(stats["CDS_length"])/1.0/len(stats["CDS_length"])
    else:
        mean_cds_num = mean_cds_length = 0

    if stats["intron_num"]:
        mean_intron_num = sum(stats["intron_num"])/1.0/len(stats["intron_num"])
        mean_intron_length =  sum(stats["intron_length"])/1.0/len(stats["intron_length"])
    else:
        mean_intron_num = mean_intron_length = 0
    if stats["exon_num"]:
       mean_exon_num = sum(stats["exon_num"])/1.0/len(stats["exon_num"])
       mean_exon_length =  sum(stats["exon_length"])/1.0/len(stats["exon_length"])
    else:
       mean_exon_num = mean_exon_length = 0

    return "{0:<10,}\t{1:<10,}\t{2:<10,}\t{3:<10,}\t{4:<10,}\t{5:<10,}\t{6:<10,}".format(stats["gene_num"], float("%.2f" % mean_gene_length), float("%.2f" % mean_cds_length), float("%.2f" % mean_exon_num), float("%.2f" % mean_exon_length), float("%.2f" % mean_intron_num), float("%.2f" % mean_intron_length))


def plot_image(labels, data_dict, outdir):
    plot_cfg = {
        "gene_length": {"window": 200, "x_max": 20000, "x_min": 0, "title": "Distribution of gene length", "x_label": "Gene length (bp)", "y_label": "Percentage (%)"},
        "CDS_length": {"window": 200, "x_max": 8000, "x_min": 0, "title": "Distribution of CDS length", "x_label": "CDS length per gene (bp)", "y_label": "Percentage (%)"},
        "exon_length": {"window": 50, "x_max": 2000, "x_min": 0, "title": "Distribution of exon length", "x_label": "Exon length (bp)", "y_label": "Number"},
        "intron_length": {"window": 50, "x_max": 2000, "x_min": 0, "title": "Distribution of intron length", "x_label": "Intron length (bp)", "y_label": "Percentage (%)"},
        "exon_num": {"window": 1, "x_max": 30, "x_min": 0, "title": "Distribution of exon number", "x_label": "Exon number per gene", "y_label": "Percentage (%)"},
        "intron_num": {"window": 1, "x_max": 30, "x_min": 0, "title": "Distribution of intron number", "x_label": "Intron number per gene", "y_label": "Percentage (%)"}
    }

    for p in plot_cfg:

        data = []
        for label in labels:
            d = data_dict[label]
            with open (os.path.join(outdir, "%s.%s.len" % (label, p)), "w") as out:
                 out.write("\n".join(map(str, sorted(d[p], reverse=True))))
            data.append(d[p])
        print("plot %s %s" % (p, labels))
        plotDistribution(data, label=labels, prefix=outdir+"/"+p, **plot_cfg[p])


def read_cfg(fn):
    r = []
    with open(fn) as fh:
         for line in fh.readlines():
             line = line.strip()
             if len(line) <=1:
                 continue

             lines = line.split()
             r.append(lines)

    return r


def main(args):
    fn = args.input
    outdir = args.outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    out = open("%s/gff.stat" % outdir, "w")
    if args.gff:
        stat_dict = gff_stat(fn, outdir)
        name = os.path.basename(fn)
        data = {name: stat_dict}
        names = [name]
        #plot_image(names, data, outdir)
        out.write("""
Number    \tAverage   \tAverage   \tAverage   \tAverage   \tAverage   \tAverage
of        \tgene      \tCDS       \texons     \texon      \tintrons   \tintro
gene      \tlength(bp)\tlength(bp)\tper gene  \tlength(bp)\tper gene  \tlength(bp)
""")
        out.write(stat2str(stat_dict))
        plot_image(names, data, outdir)
    if args.config:
        cfg = read_cfg(fn)
        out.write("""
Source         \tNumber    \tAverage   \tAverage   \tAverage   \tAverage   \tAverage   \tAverage
               \tof        \tgene      \tCDS       \texons     \texon      \tintrons   \tintro
               \tgene      \tlength(bp)\tlength(bp)\tper gene  \tlength(bp)\tper gene  \tlength(bp)
""")
        data_dict = {}
        names = []
        for s in cfg:
            name = s[0]
            fn = s[1]
            print("Statstic for species %s" % name)
            stat_dict = gff_stat(fn, outdir)
            out.write("%-15s\t%s\n" % (name, stat2str(stat_dict)))
            data_dict[name] = stat_dict
            names.append(name)
        plot_image(names, data_dict, outdir)
    out.close()


if __name__ == "__main__":
    opt = argparse.ArgumentParser(formatter_class=HelpFormatter,description="""
description:
    Statstic of gene, CDS, exon, intron in gff3

    config example:
    E.coli      /path/ecol.gff3

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    group = opt.add_mutually_exclusive_group(required=True)
    group.add_argument("-g", "--gff", action="store_true",
                    help="gff3 format input")
    group.add_argument("-c", "--config", action="store_true",
                    help="config format 'species_name  gff_path'" )
    group.add_argument("-o", "--outdir", default="gffStat.out",
                    help="output dir" )
    opt.add_argument("input",
                    help="the input filename" )

    args = opt.parse_args()
    main(args)


