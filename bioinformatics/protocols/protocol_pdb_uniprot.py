# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import lxml.etree as ET
import os
import sys
import urllib.request

from pwem.protocols import EMProtocol
import pyworkflow.object as pwobj
from pyworkflow.protocol.params import (PointerParam)
from bioinformatics.objects import DatabaseID, SetOfDatabaseID

class ProtBioinformaticsPDBUniprot(EMProtocol):
    """Query PDB for Uniprot sequences related to these proteins"""

    """
        Retrieve UniProt sequence annotations associated with Protein
        Data Bank (PDB) structures. The protocol queries the RCSB PDB
        mapping services, identifies corresponding UniProt entries, and
        generates a structured collection of protein sequence references
        for downstream structural biology and bioinformatics workflows.

        AI Generated:

        PDB to UniProt (ProtBioinformaticsPDBUniprot) — User Manual
            Overview

            The PDB to UniProt protocol is designed to retrieve UniProt
            sequence annotations associated with experimentally resolved
            Protein Data Bank structures. Its primary objective is to
            establish a direct relationship between macromolecular
            structural entries and their corresponding protein sequence
            records available in the UniProt database.

            In practical structural biology and bioinformatics workflows,
            this protocol is particularly useful for integrating
            experimentally determined structural data with curated protein
            sequence annotations, functional descriptions, and biological
            classifications. By automatically querying the RCSB PDB
            mapping services, the protocol simplifies the identification
            and organization of UniProt references linked to specific PDB
            structures.

            Inputs and General Workflow

            The protocol requires a collection of PDB identifiers provided
            as a SetOfDatabaseID object. During execution, the protocol
            iterates through the supplied structural entries and queries
            the RCSB PDB UniProt mapping service for each structure.

            For each input entry, the protocol generates a local XML
            response containing alignment and mapping information between
            the queried PDB structure and UniProt sequence records. If
            chain identifiers are available in the input object, the
            protocol incorporates chain-specific queries to improve the
            precision of sequence mapping in multi-chain assemblies.

            Retrieved XML responses are parsed to identify UniProt
            accession identifiers associated with the structural entry.
            The protocol then enriches the original database objects with
            UniProt accession information and direct links to the
            corresponding UniProt repository entries.

            UniProt Mapping and Biological Annotation

            For each successfully identified UniProt record, the protocol
            stores the UniProt accession identifier together with a direct
            URL pointing to the associated UniProt entry. These mappings
            establish a biologically meaningful bridge between structural
            and sequence-based information resources.

            From a biological perspective, this integration is especially
            valuable because UniProt records contain curated annotations
            including protein function, sequence conservation, domain
            organization, post-translational modifications, taxonomic
            information, and known biological variants.

            The protocol is therefore highly useful when preparing
            datasets for comparative sequence analysis, structural
            annotation, evolutionary studies, integrative modeling, or
            functional interpretation of experimentally resolved
            macromolecular structures.

            Error Handling and Data Availability

            If no UniProt mapping is identified for a given PDB entry, the
            protocol automatically assigns placeholder values indicating
            that UniProt information is not available. This mechanism
            preserves metadata consistency across all generated entries
            and prevents incomplete object definitions during downstream
            processing.

            XML parsing and remote retrieval failures are handled
            gracefully, allowing the protocol to continue processing
            additional structures even when isolated queries cannot be
            successfully resolved.

            Outputs and Interpretation

            The protocol produces a SetOfDatabaseID object containing the
            original structural entries enriched with UniProt-related
            metadata. Each output entry preserves the associated PDB
            information while incorporating the mapped UniProt accession
            and corresponding external reference link.

            These outputs can subsequently be integrated into structural
            annotation workflows, protein classification analyses,
            comparative bioinformatics studies, sequence retrieval
            pipelines, or downstream structural and functional modeling
            applications.

            Practical Recommendations

            In routine structural biology workflows, this protocol is
            particularly useful when linking experimentally resolved
            structures with curated sequence repositories, organizing
            datasets according to protein families, or integrating
            structural and functional annotations into unified biological
            analyses.

            Since the protocol depends on online access to the RCSB PDB
            infrastructure, stable internet connectivity is important for
            reliable execution. For large collections of structures,
            preserving locally downloaded XML responses may improve
            reproducibility and reduce repeated network access.

            Final Perspective

            The PDB to UniProt protocol provides an efficient bridge
            between experimentally determined structural data and curated
            protein sequence annotations. By automatically retrieving and
            organizing UniProt mappings associated with PDB structures,
            the protocol facilitates biologically meaningful integration
            of structural and sequence information within structural
            biology and bioinformatics workflows.
        """
    _label = 'pdb -> uniprot'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputListID', PointerParam, pointerClass="SetOfDatabaseID",
                       label='List of PDB Ids:', allowsNull=False,
                       help="List of atomic structures for the query")

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('searchStep')

    def searchStep(self):
        outputDatabaseID = SetOfDatabaseID().create(path=self._getPath())
        for item in self.inputListID.get():
            newItem = DatabaseID()
            newItem.copy(item)
            newItem._uniprotId = pwobj.String("Not available")
            newItem._uniprotLink = pwobj.String("Not available")

            pdbId = item._pdbId.get()
            print("Processing %s"%pdbId)

            urlId = "https://www.rcsb.org/pdb/rest/das/pdb_uniprot_mapping/alignment?query=%s" % pdbId
            if hasattr(item,"_chain"):
                urlId+="."+item._chain.get().upper()

            fnXml=self._getExtraPath("%s.xml"%pdbId)
            if not os.path.exists(fnXml):
                print("Fetching uniprot: %s"%urlId)
                for i in range(3):
                    try:
                        urllib.request.urlretrieve(urlId,fnXml)
                        break
                    except: # The library raises an exception when the web is not found
                        pass
            if os.path.exists(fnXml):
                try:
                    tree = ET.parse(fnXml)
                    # print(ET.tostring(tree, pretty_print=True))

                    uniprotId = None
                    for child in tree.getroot().iter():
                        if child.tag.endswith("alignObject"):
                            if child.attrib['dbSource']=="UniProt":
                                uniprotId=child.attrib['dbAccessionId']
                                break
                    if uniprotId:
                        newItem._uniprotId = pwobj.String(uniprotId)
                        newItem._uniprotLink = pwobj.String("https://www.uniprot.org/uniprot/%s"%uniprotId)
                except:
                    print("    Cannot parse the Uniprot XML: %s"%fnXml)

                outputDatabaseID.append(newItem)

        self._defineOutputs(outputUniprot=outputDatabaseID)
        self._defineSourceRelation(self.inputListID, outputDatabaseID)
