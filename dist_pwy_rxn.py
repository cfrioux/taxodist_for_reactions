#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 09:05:00 2019

@author: cfrioux
"""

from padmet.classes import PadmetSpec
import time
import json
import argparse
import sys
from json import JSONEncoder
from taxon import Taxon, TaxonNotFound

class MyEncoder(JSONEncoder):
        def default(self, o):
            return o.__dict__

def pathway_info(metacycdb):
    """Retrieve pathways from Metacyc.
    
    Reads a metacyc padmet file and return a set of taxa as well as the
    association between reactions and taxa (via the belonging of reactions
    in pathways).

    Args:
        metacycdb (str): Metacyc padmet file
    
    Returns:
        set,dict: taxonomy info of rxn, taxons occuring in all pathways
    """
    padmet = PadmetSpec(metacycdb)

    rxn_in_pwy = {}
    all_taxa = set()

    for rlt in padmet.getAllRelation():
        if rlt.type == "is_in_pathway" and "TAXONOMIC-RANGE" in padmet.dicOfNode[rlt.id_out].misc.keys():
            rxn_id, pwy_id = rlt.id_in, rlt.id_out
            taxons_set = set(padmet.dicOfNode[pwy_id].misc["TAXONOMIC-RANGE"])
            all_taxa.update(taxons_set)
            try:
                rxn_in_pwy[rxn_id].update(taxons_set)
            except KeyError:
                rxn_in_pwy[rxn_id] = taxons_set
    return all_taxa, rxn_in_pwy


def print_taxa_of_rxn(rxn_file, dict_rxn):
    """Print the taxa associated to reactions in a file.
    
    Args:
        rxn_file (str): rxn file, one rxn per line
        dict_rxn (dict): taxonomy info of rxn
    """
    with open(rxn_file,"r") as f:
        reactions = f.read().splitlines()
    for reaction in reactions:
        try:
            print(dict_rxn[reaction])
        except KeyError:
            print("no associated pathway")
    return

def turn_taxa_into_taxonomic_info(set_of_taxa, email=None):
    """Retrieve taxonomic information for all taxa in a set with NCBI api.
    
    Args:
        set_of_taxa (set): set of taxa

    Returns:
        dict: dictionary of taxon objects for all IDs of the input set
    """
    tax_info = {}
    for elem in set_of_taxa:
        # avoid more than 3 requests per seconds so wait 0.35s between each call to NCBI
        tax_number = elem.strip('TAX-')
        try:
            taxon = Taxon(taxid=tax_number, useremail=email)
            tax_info[elem] = taxon
        except TaxonNotFound:
            print("Taxonomic ID " + tax_number + " not found")
            pass
        time.sleep(.350)
    return(tax_info)

def write_json(data, jsonfile):
    """Write a json file.
    
    Args:
        data (dict): dictionary key = taxon str, value = Taxon
        jsonfile (str): file name for json
    """
    outjson = open(jsonfile, 'w')
    outjson.write(json.dumps(data, cls=MyEncoder))
    outjson.close()

def read_json(jsonfile):
    """Load dict from json file.
    
    Args:
        jsonfile (str): json file
    
    Returns:
        dict: data
    """
    with open(jsonfile, 'r') as f:
        data = json.load(f)
    data_obj = {}
    for elem in data:
        data_obj[elem] = Taxon(taxid=data[elem]["taxid"],
                scientific_name=data[elem]["scientific_name"],
                lineage_taxa_name=data[elem]["lineage_taxa_name"],
                lineage_taxa_id=data[elem]["lineage_taxa_name"],
                parent_taxid=data[elem]["parent_taxid"])
    return data_obj

def distance_org_taxa(organism, taxaset, taxodic, email=None, penalty_price = 20):
    """Calculate distance between given organism and any of the taxa in a set.
    
    Args:
        organism (str): name of the organism
        taxaset (set): taxa
        taxodic (dict): dictionary of taxo information
        penalty_price (int, optional): Defaults to 20. penalty for changing branch in taxo tree
    
    Returns:
        dict: taxo distance of organism to every taxon
    """
    target_orga = Taxon(scientific_name=organism, useremail=email)
    taxdist_dic = {}
    for taxa in taxaset:
        try:
            taxdist_dic[taxa] = target_orga.get_distance_between_two_taxa(taxodic[taxa])
        except KeyError:
            print(taxa + " is not in the taxo data. Ignoring it. Either it \
was not in the NCBI datasbae, or the json file is outdated wrt \
the padmet reference database."                                                                                                                                                                                                                                                                                                                                                       )
    return taxdist_dic

def associate_rxn_dist(rxndata, taxdistdata):
    """Associate a distance to each rxn
    
    Args:
        rxndata (dict): association of taxa for every rxn
        taxdistdata (dict): distance of orga to every taxa
    
    Returns:
        dict: distance for all rxn 
    """
    all_rxn_dist = {}
    for rxn in rxndata:
        try:
            rxn_dist = [taxdistdata[txn] for txn in rxndata[rxn]]
        except KeyError: #there might be taxa in Metacyc db that were not found in NCBI
            rxn_dist = []
            for txn in rxndata[rxn]:
                try:
                    rxn_dist.append(taxdistdata[txn])
                except KeyError:
                    pass
        try:
            all_rxn_dist[rxn] = min(rxn_dist)
        except ValueError: #no real NCBI taxa was found for this rxn
            all_rxn_dist[rxn] = 1000
    return all_rxn_dist


if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage="dist_pwy_rxn.py --email [--padmet] [--fromjson] [--tojson] [--orga]")
    parser.add_argument("--email", dest="email", help = "email is required to call NCBI API", required=True)
    parser.add_argument("--orga", dest="orga", help = "organism of interest", required=False)
    parser.add_argument("--fromjson", dest="fromjson", help = "json input to retrieve data", required=False)
    parser.add_argument("--tojson", dest="tojson", help = "json output to store data", required=False)
    parser.add_argument("--padmet", dest="padmet", help = "padmet file with raw data on taxonomy and pathways", required=False)

    args = parser.parse_args()
    if len(sys.argv[1:]) == 1:
        parser.print_help()
        parser.exit()
            
    email = args.email
    
    if args.fromjson and args.orga and args.padmet:
        taxoinfo = read_json(args.fromjson)
        metacyc_taxa, reactions_asso_taxa = pathway_info(args.padmet)
        taxodist = distance_org_taxa(args.orga, metacyc_taxa, taxoinfo, email)
        rxn_asso_dist = associate_rxn_dist(reactions_asso_taxa, taxodist)
        print(len(rxn_asso_dist))
        
    elif args.padmet and args.tojson and args.orga:
        metacyc_taxa, reactions_asso_taxa = pathway_info(args.padmet)
        taxoinfo = turn_taxa_into_taxonomic_info(metacyc_taxa, email)
        write_json(taxoinfo, args.tojson)
        taxodist = distance_org_taxa(args.orga, metacyc_taxa, taxoinfo, email)
        rxn_asso_dist = associate_rxn_dist(reactions_asso_taxa, taxodist)
        print(len(rxn_asso_dist))
        
    elif args.padmet and args.orga:
        metacyc_taxa, reactions_asso_taxa = pathway_info(args.padmet)
        taxoinfo = turn_taxa_into_taxonomic_info(metacyc_taxa, email)
        taxodist = distance_org_taxa(args.orga, metacyc_taxa, taxoinfo, email)
        rxn_asso_dist = associate_rxn_dist(reactions_asso_taxa, taxodist)
        print(len(rxn_asso_dist))
        
    elif args.padmet and args.tojson:
        metacyc_taxa, reactions_asso_taxa = pathway_info(args.padmet)
        write_json(turn_taxa_into_taxonomic_info(metacyc_taxa, email), args.tojson)
    
    else:
        parser.print_help()
        parser.exit()
        

