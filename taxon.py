from Bio import Entrez
import logging
from urllib.error import HTTPError

logger = logging.getLogger(__name__)

class TaxonNotFound(Exception):
    """Define a taxon not found.
    """

    def __init__(self, taxdata):
        """Construct an exception in case taxonomic data is not found in NCBI.
        
        Args:
            taxtaxdataid (str): taxonomic data (ID or scientific name)        
        """
        super(TaxonNotFound,self).__init__()
        logger.error(f"Taxonomic information for {taxdata} not found in NCBI")


class TaxonInvalidInput(Exception):
    """Define an exception for invalid input.
    """

    def __init__(self):
        """Construct an exception in case input is invalid.
        """
        super(TaxonInvalidInput, self).__init__()
        logger.error("Invalid input to create a taxon")


class Taxon():
    """Define a taxon object from NCBI API

    >>> Taxon("2",useremail="nobody@example.com").scientific_name
    'Bacteria'
    >>> Taxon("blop",useremail="nobody@example.com")
    Traceback (most recent call last):
        ...
    TaxonNotFound
    >>> Taxon(scientific_name="Ectocarpus siliculosus",useremail="nobody@example.com").taxid
    '2880'
    """

    def __init__(self, taxid: str=None, scientific_name: str=None, 
                 lineage_taxa_name: list=None, lineage_taxa_id: list=None, 
                 parent_taxid: str=None, useremail: str=None):
        """Construct a taxon with NCBI API.
    
        Args:
            taxid (str, optional): NCBI taxonomy ID
            scientific_name (str, optional): NCBI scientific name
            lineage_taxa_name (list, optional): list of all parents names
            lineage_taxa_id (list, optional): list of all parents IDs
            parent_taxid (str, optional): ID of parent taxon
            useremail (str, optional): email to call NCBI API

        Raises:
            TaxonNotFound
        """
        #first way to build object: all info is given, eg from reading a json
        if taxid and scientific_name and lineage_taxa_name and parent_taxid and lineage_taxa_id:
            self.taxid = taxid
            self.scientific_name = scientific_name
            self.lineage_taxa_name = lineage_taxa_name
            self.lineage_taxa_id = lineage_taxa_id
            self.parent_taxid = parent_taxid
            return
        else:
            if useremail:
                Entrez.email = useremail  # Always tell NCBI who you are
            else:
                logger.warning("It is strongly advised to use an email to connect to the NCBI API")
            if not taxid and not scientific_name:
                raise TaxonInvalidInput
            #second way to build object: having name and calling API to get ID, and continue
            elif not taxid:
                self.taxid = self.get_taxid_from_name(scientific_name)
            #third and last way to build object: having ID and calling API
            else:
                self.taxid = taxid
            try:
                handle = Entrez.efetch(db="taxonomy", id=self.taxid, rettype="xml")
                record = Entrez.read(handle)
                entry = record[0]
            except (TypeError, IndexError, HTTPError):
                raise TaxonNotFound(self.taxid)
            self.scientific_name = entry.get('ScientificName')
            self.lineage_taxa_name = entry.get('Lineage').split('; ')
            self.lineage_taxa_name.append(self.scientific_name)
            self.parent_taxid = entry['ParentTaxId']
            try:
                self.lineage_taxa_id = [x['TaxId'] for x in entry['LineageEx']]
                self.lineage_taxa_id.append(self.taxid)  #add taxid at the end of lineage
            except KeyError: # cellular organisms has no LineageEx
                self.lineage_taxa_id = [self.taxid]
            #calculate distance to the root organism
            if self.taxid==1:
                self.dist_to_root = 0
            else:
                self.dist_to_root = len(self.lineage_taxa_id)
            
    def get_taxid_from_name(self, scientific_name: str):
        """Get taxonomic id from scientific name
        
        Args:
            scientific_name (str): scientific name
        """
        handle = Entrez.esearch(
            db="taxonomy", term=scientific_name, rettype="xml")
        record = Entrez.read(handle)
        try:
            return record["IdList"][0]
        except (TypeError, IndexError):
            raise TaxonNotFound(scientific_name)

    def get_distance_between_two_taxa(self, anothertaxon: 'Taxon', penalty_price: int=20):
        """Get distance between a Taxon object and another.
        
        Args:
            anothertaxon (Taxon): a Taxon
            penalty_price (int, optional): Defaults to 20. price of penalty to branch the tree
        
        Returns:
            int: distance between the two taxa

        >>> tax1 = Taxon(taxid="2880", useremail="nobody@example.com")
        >>> tax2 = Taxon(taxid="2", useremail="nobody@example.com")
        >>> tax1.get_distance_between_two_taxa(tax2)
        28
        """
        if anothertaxon.taxid == 1:
            return len(self.lineage_taxa_id)
        else:
            common_tax_id = []
            down_distance = 0
            up_distance = 0
            i = 0
            if len(self.lineage_taxa_id) <= len(anothertaxon.lineage_taxa_id):
                while i < len(self.lineage_taxa_id):
                    if self.lineage_taxa_id[i] == anothertaxon.lineage_taxa_id[
                            i]:
                        common_tax_id.append(self.lineage_taxa_id[i])
                        i = i+1
                    else:
                        break
            else:
                while i < len(anothertaxon.lineage_taxa_id):
                    if self.lineage_taxa_id[i] == anothertaxon.lineage_taxa_id[
                            i]:
                        common_tax_id.append(self.lineage_taxa_id[i])
                        i = i+1
                    else:
                        break
            #print(len(common_tax_id))
            # below should be used for a unit test
            #if common_tax_id == []:
            #    print(id_txn)
            up_distance = len(self.lineage_taxa_id) - len(common_tax_id)
            # distance from organism to a node higer in the tree :
            # the youngest common ancestor with the expected taxonomic
            # range or root if there is no common ancestor
            if up_distance != 0:
                down_distance = len(
                    anothertaxon.lineage_taxa_id) - len(common_tax_id)
                #distance from top of the tree (around root or so)
                #to bottom (taxon)
            else:
                down_distance = 0
            return up_distance + down_distance * penalty_price


if __name__ == "__main__":
    import doctest
    doctest.testmod()
