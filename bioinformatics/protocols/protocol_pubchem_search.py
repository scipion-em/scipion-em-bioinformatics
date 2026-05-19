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

import contextlib
import lxml.etree as ET
import os
import urllib.request

import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, BooleanParam, EnumParam

class ProtBioinformaticsPubChemSearch(EMProtocol):
    """Add the best batching entry from Pubchem https://pubchem.ncbi.nlm.nih.gov/"""

    """
        Retrieve and annotate small molecules with PubChem information
        using chemical identifiers derived from molecular structure files.
        The protocol queries the PubChem database, identifies matching
        compounds through SMILES representations, and enriches molecular
        entries with standardized PubChem annotations for downstream
        cheminformatics and bioinformatics workflows.

        AI Generated:

        PubChem Search (ProtBioinformaticsPubChemSearch) — User Manual
            Overview

            The PubChem Search protocol is designed to associate small
            molecules with their corresponding PubChem compound entries.
            Its primary purpose is to enrich molecular datasets with
            standardized chemical annotations obtained directly from the
            PubChem repository.

            In practical cheminformatics and structural biology workflows,
            this protocol is especially useful for identifying compounds,
            validating molecular annotations, and connecting locally stored
            molecules with publicly available chemical information. By
            integrating PubChem metadata into Scipion workflows, the
            protocol facilitates interoperability between molecular
            datasets and external bioinformatics resources.

            Inputs and General Workflow

            The protocol requires a SetOfSmallMolecules object containing
            molecular structure files. During execution, each molecule is
            analyzed individually and its SMILES representation is used to
            query the PubChem REST API.

            When a matching compound is identified, the protocol retrieves
            the associated PubChem Compound ID (CID) together with the
            compound name and stores this information locally for reuse in
            subsequent executions. This caching mechanism reduces repeated
            network requests and improves workflow reproducibility.

            Compound Annotation and Metadata Retrieval

            For each successfully matched molecule, the protocol creates an
            updated molecular entry enriched with PubChem annotations. The
            retrieved metadata includes the compound name and a direct URL
            linking to the corresponding PubChem entry.

            From a biological and chemical perspective, this functionality
            is valuable for compound identification, ligand annotation,
            molecular cataloging, and integration with external chemical
            databases. The protocol simplifies the connection between local
            molecular structures and globally recognized chemical
            repositories.

            Outputs and Interpretation

            The protocol produces an updated set of small molecules where
            each entry contains additional PubChem annotations when a valid
            match is found. Molecules without successful matches are still
            preserved in the output dataset, ensuring dataset consistency.

            These outputs can subsequently be used in molecular screening,
            ligand analysis, annotation workflows, or integration with
            external cheminformatics pipelines requiring PubChem
            identifiers and references.

            Final Perspective

            The PubChem Search protocol provides an efficient mechanism for
            connecting molecular structure datasets with standardized
            chemical information available in PubChem. By automatically
            retrieving compound identifiers and annotations, the protocol
            supports biologically meaningful molecular interpretation and
            improves interoperability across bioinformatics and
            cheminformatics workflows.
        """
    _label = 'PubChem'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Set to filter:', allowsNull=False)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('operateStep')

    def operateStep(self):
        outputSet = self.inputSet.get().create(self._getPath())
        for oldEntry in self.inputSet.get():
            fnSmall = oldEntry.smallMoleculeFile.get()
            fnBase = os.path.splitext(os.path.split(fnSmall)[1])[0]
            fnName = self._getExtraPath(fnBase+".txt")
            cid = None
            if not os.path.exists(fnName):
                if fnSmall.endswith('.smi'):
                    fhSmile = open(fnSmall)
                    smile = fhSmile.readlines()[0].split()[0].strip() # Only first line
                    fhSmile.close()
                    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/%s/cids/TXT"%smile
                    print(url)

                pubChemName = ""
                try:
                    with contextlib.closing(urllib.request.urlopen(url)) as fp:
                        mybytes = fp.read()
                        cid = mybytes.decode("utf8").split()[0]
                        fp.close()

                    if cid!="0":
                        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%s/XML/?response_type=save&response_basename=compound_CID_%s"%(cid,cid)
                        print(url)
                        fnXml = self._getTmpPath("compound.xml")
                        urllib.request.urlretrieve(url, fnXml)
                        if os.path.exists(fnXml):
                            tree = ET.parse(fnXml)

                        pubChemName = ""
                        for child in tree.getroot().iter():
                            if "RecordTitle" in child.tag:
                                pubChemName = child.text
                                print(pubChemName)
                                fh = open(fnName,'w')
                                fh.write(pubChemName+" ;; %s"%cid)
                                fh.close()
                                break
                except Exception as e:
                    print(e)
                    print("  Could not be retrieved")
            else:
                fh = open(fnName)
                tokens = fh.readline().split(';;')
                pubChemName = tokens[0].strip()
                cid = tokens[1].strip()
                fh.close()

            newEntry = self.inputSet.get().ITEM_TYPE()
            newEntry.copy(oldEntry)
            newEntry.pubChemName = pwobj.String(pubChemName)
            if cid is not None and cid!="0":
                url = "https://pubchem.ncbi.nlm.nih.gov/compound/%s"%cid
            else:
                url = ""
            newEntry.pubChemURL = pwobj.String(url)
            outputSet.append(newEntry)

        if len(outputSet)>0:
            self._defineOutputs(output=outputSet)
            self._defineSourceRelation(self.inputSet, outputSet)

    def _validate(self):
        errors = []
        firstItem = self.inputSet.get().getFirstItem()
        if not hasattr(firstItem,"smallMoleculeFile"):
            errors.append("The input set does not contain small molecules")
        return errors
