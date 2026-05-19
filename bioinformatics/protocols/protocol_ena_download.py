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
import urllib.request

from pwem.protocols import EMProtocol
from pwem.convert.sequence import sequenceLength
import pyworkflow.object as pwobj
from pyworkflow.protocol.params import PointerParam
from bioinformatics.objects import DatabaseID, SetOfDatabaseID, NucleotideSequenceFile

class ProtBioinformaticsEnaDownload(EMProtocol):
    """Download the Fasta files of a set of enaId's"""

    """
        Download nucleotide FASTA sequences from the ENA database using
        ENA identifiers associated with input database entries.

        AI Generated:

        ENA Sequence Download (ProtBioinformaticsEnaDownload) — User Manual
            Overview

            The ENA Sequence Download protocol retrieves nucleotide
            sequences from the European Nucleotide Archive (ENA) using
            ENA accession identifiers stored in a set of database entries.
            Its main purpose is to automate sequence retrieval and
            integrate nucleotide data into structural biology and
            bioinformatics workflows.

            In many biological studies, nucleotide sequences are required
            for sequence comparison, annotation, evolutionary analysis,
            or integrative modeling. This protocol simplifies the process
            of obtaining ENA sequence data while maintaining the
            association between downloaded sequences and their original
            biological identifiers.

            Inputs and Workflow

            The protocol requires a set of database identifiers containing
            valid ENA accession IDs. During execution, each identifier is
            queried against the ENA repository through its public web
            service interface.

            For every valid entry, the corresponding FASTA sequence is
            downloaded and stored locally. The protocol also calculates
            sequence length information and associates it with each
            retrieved entry. All successfully downloaded sequences are
            additionally merged into a combined FASTA file for easier
            downstream analysis.

            Outputs and Interpretation

            After execution, the protocol generates an updated set of
            database entries enriched with sequence file references and
            sequence metadata. A combined FASTA file containing all
            retrieved nucleotide sequences is also produced.

            These outputs can be directly used in downstream
            bioinformatics applications such as sequence alignment,
            annotation, phylogenetic analysis, or integrative structural
            biology workflows.

            Final Perspective

            The protocol provides a practical and automated solution for
            integrating ENA nucleotide sequence retrieval into biological
            and bioinformatics pipelines, facilitating reproducible and
            efficient sequence-based analyses.
        """
    _label = 'ena download'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputListID', PointerParam, pointerClass="SetOfDatabaseID",
                       label='List of Uniprot Ids:', allowsNull=False,
                       help="List of atomic structures for the query")

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('searchStep')

    def searchStep(self):
        outputDatabaseID = SetOfDatabaseID().create(path=self._getPath())
        fnList = []
        for item in self.inputListID.get():
            newItem = DatabaseID()
            newItem.copy(item)
            newItem._enaFile = pwobj.String("Not available")
            newItem._enaSeqLength = pwobj.Integer(-1)

            enaId = item._enaId.get()
            print("Processing %s"%enaId)

            urlId = "https://www.ebi.ac.uk/ena/data/view/%s&display=fasta" % enaId

            fnFasta=self._getExtraPath("%s.fasta"%enaId)
            if not os.path.exists(fnFasta):
                print("Fetching ena: %s"%urlId)
                for i in range(3):
                    try:
                        urllib.request.urlretrieve(urlId,fnFasta)
                        if not fnFasta in fnList:
                            fnList.append(fnFasta)
                        break
                    except: # The library raises an exception when the web is not found
                        pass
            if os.path.exists(fnFasta):
                newItem._enaFile = pwobj.String(fnFasta)
                newItem._enaSeqLength = pwobj.Integer(sequenceLength(fnFasta))

            outputDatabaseID.append(newItem)

        fnAll = self._getPath("sequences.fasta")
        with open(fnAll, 'w') as outfile:
            for fname in fnList:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
                    outfile.write('\n\n')
        seqFile=NucleotideSequenceFile()
        seqFile.setFileName(fnAll)

        self._defineOutputs(outputUniprot=outputDatabaseID)
        self._defineSourceRelation(self.inputListID, outputDatabaseID)
        self._defineOutputs(outputSequence=seqFile)
        self._defineSourceRelation(self.inputListID, seqFile)

    def _validate(self):
        errors=[]
        if not hasattr(self.inputListID.get().getFirstItem(),"_enaId"):
            errors.append("The set does not have an _enaId")
        return errors