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

import os
import sys
import urllib.request

from pwem.protocols import EMProtocol
import pyworkflow.object as pwobj
from pyworkflow.protocol.params import (PointerParam)
from bioinformatics.objects import DatabaseID, SetOfDatabaseID

class ProtBioinformaticsSmallMoleculesPDB(EMProtocol):
    """Extract a set of PDB Ids from the interaction list of a set of small molecules"""

    """
        Extract and organize Protein Data Bank (PDB) identifiers associated
        with a collection of small molecules. The protocol analyzes stored
        interaction metadata and generates a structured set of PDB entries
        linked to the corresponding interacting compounds for downstream
        structural biology and bioinformatics workflows.

        AI Generated:

        Small Molecules to PDB (ProtBioinformaticsSmallMoleculesPDB) — User Manual
            Overview

            The Small Molecules to PDB protocol is designed to identify
            Protein Data Bank structures associated with a collection of
            small molecules. Its primary objective is to transform
            ligand-centered interaction information into a structure-
            centered dataset that can be used for structural analysis,
            interaction studies, or database exploration workflows.

            In practical bioinformatics applications, this protocol is
            useful when studying ligand-binding proteins, identifying
            structural targets associated with specific compounds, or
            exploring interaction relationships between chemical entities
            and experimentally solved macromolecular structures.

            Inputs and General Workflow

            The protocol requires a set of database identifiers
            representing small molecules enriched with interaction
            annotations. During execution, the protocol scans the
            interaction metadata stored within each molecule and extracts
            the associated PDB identifiers.

            The retrieved PDB entries are consolidated into a unique
            structure collection while preserving information about the
            interacting compounds linked to each structure. This allows
            the protocol to generate a biologically meaningful mapping
            between ligands and structural entries.

            Interaction Mapping and Metadata Organization

            For each detected PDB identifier, the protocol creates a
            database entry containing the corresponding structure code,
            direct access links to the RCSB PDB repository, and the list
            of interacting chemical identifiers associated with that
            structure.

            When available, the protocol also supports PDB chemical
            identifiers, enabling compatibility with workflows that use
            ligand annotations derived directly from PDB chemical
            component definitions.

            Outputs and Interpretation

            The protocol produces a SetOfDatabaseID object containing the
            extracted PDB structures and their associated ligand
            interaction metadata. Each output entry represents a unique
            structural record connected to one or more interacting small
            molecules.

            These outputs can subsequently be integrated into structural
            annotation workflows, interaction network analyses, ligand-
            binding studies, or downstream bioinformatics pipelines.

            Final Perspective

            The Small Molecules to PDB protocol provides an efficient
            bridge between ligand-centered datasets and structural
            repositories. By automatically organizing interaction-derived
            PDB identifiers into a structured collection, the protocol
            facilitates biologically meaningful exploration of protein–
            ligand relationships within structural biology workflows.
        """
    _label = 'small mols -> pdb'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputListID', PointerParam, pointerClass="SetOfDatabaseID",
                       label='List of small molecules:', allowsNull=False,
                       help="List of small molecules to extract the list of PDB IDs")

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('extractStep')

    def extractStep(self):
        listIds={}
        smallMolID = "none"
        for item in self.inputListID.get():
            if hasattr(item,"_iteractsWithPDBId"):
                tokens = item._iteractsWithPDBId.get().split(";")
                for token in tokens:
                    pdbId = token.strip()
                    if not pdbId in listIds:
                        listIds[pdbId]=[]
                    chemId=item.getDbId()
                    if hasattr(item,"_PDBChemId"):
                        chemId=item._PDBChemId.get()
                        smallMolID = "pdbchem"
                    listIds[pdbId].append(chemId)

        outputDatabaseID = SetOfDatabaseID().create(path=self._getPath(),suffix='PDBs')
        for pdbId in listIds:
            pdb = DatabaseID()
            pdb.setDatabase("pdb")
            pdb.setDbId(pdbId)
            pdb._pdbId = pwobj.String(pdbId)
            pdb._PDBLink = pwobj.String("https://www.rcsb.org/structure/%s" % pdbId)
            aux = " ; ".join(listIds[pdbId])
            if smallMolID == "none":
                pdb._interactsWithChemId = pwobj.String(aux)
            elif smallMolID == "pdbchem":
                pdb._interactsWithPDBChemId = pwobj.String(aux)
            outputDatabaseID.append(pdb)
        self._defineOutputs(outputPDBs=outputDatabaseID)
        self._defineSourceRelation(self.inputListID, outputDatabaseID)
