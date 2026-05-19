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
from pwem.convert.sequence import sequenceLength
import pyworkflow.object as pwobj
from pyworkflow.protocol.params import PointerParam
from bioinformatics.objects import DatabaseID, SetOfDatabaseID, ProteinSequenceFile

class ProtBioinformaticsUniprotDownload(EMProtocol):
    """Download the Fasta files of a set of uniprotId's"""

    """
        Download and organize protein FASTA sequences associated with
        UniProt identifiers. The protocol retrieves sequence files from
        the UniProt repository and generates both individual and combined
        sequence datasets for downstream bioinformatics analyses.

        AI Generated:

        UniProt Sequence Download (ProtBioinformaticsUniprotDownload) — User Manual
            Overview

            The UniProt Sequence Download protocol is designed to retrieve
            protein FASTA sequences directly from the UniProt database using
            a collection of UniProt accession identifiers. Its primary goal
            is to automate sequence acquisition and integrate protein
            sequence information into Scipion-based bioinformatics and
            structural biology workflows.

            In practical applications, this protocol is commonly used when
            preparing datasets for sequence alignment, comparative genomics,
            homology modeling, structural prediction, or functional
            annotation analyses. By automatically downloading and organizing
            sequence information, the protocol simplifies large-scale
            sequence management tasks.

            Inputs and General Workflow

            The protocol requires a set of database identifiers containing
            UniProt accession codes. During execution, each UniProt entry is
            processed independently, and the corresponding FASTA sequence is
            retrieved from the UniProt online repository.

            Downloaded sequences are stored locally as individual FASTA
            files. The protocol also tracks sequence metadata, including the
            local file location and the total protein sequence length for
            each retrieved entry.

            Sequence Retrieval and Organization

            For each UniProt identifier, the protocol accesses the UniProt
            FASTA endpoint and attempts to download the associated protein
            sequence. To improve robustness against temporary network or
            server interruptions, multiple retrieval attempts are performed
            automatically.

            Successfully downloaded FASTA files are integrated into the
            output dataset and annotated with their corresponding sequence
            lengths. In addition, all retrieved sequences are merged into a
            single consolidated FASTA file containing the complete sequence
            collection generated during the execution.

            Outputs and Interpretation

            The protocol produces an updated SetOfDatabaseID object enriched
            with sequence-related metadata, including references to the
            downloaded FASTA files and their associated protein lengths.

            Additionally, a ProteinSequenceFile object is generated,
            containing a combined FASTA file with all retrieved protein
            sequences. This consolidated dataset can be directly used in
            downstream sequence-based workflows, including alignments,
            database searches, structural modeling, or functional analyses.

            Final Perspective

            The UniProt Sequence Download protocol provides an efficient and
            automated mechanism for integrating UniProt protein sequences
            into bioinformatics workflows. By organizing both individual and
            combined FASTA datasets together with sequence metadata, the
            protocol facilitates scalable and biologically meaningful
            sequence analysis pipelines.
        """
    _label = 'uniprot download'

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
            newItem._uniprotFile = pwobj.String("Not available")
            newItem._unitprotSeqLength = pwobj.Integer(-1)

            uniprotId = item._uniprotId.get()
            print("Processing %s"%uniprotId)

            urlId = "https://www.uniprot.org/uniprot/%s.fasta" % uniprotId

            fnFasta=self._getExtraPath("%s.fasta"%uniprotId)
            if not os.path.exists(fnFasta):
                print("Fetching uniprot: %s"%urlId)
                for i in range(3):
                    try:
                        urllib.request.urlretrieve(urlId,fnFasta)
                        if not fnFasta in fnList:
                            fnList.append(fnFasta)
                        break
                    except: # The library raises an exception when the web is not found
                        pass
            if os.path.exists(fnFasta):
                newItem._uniprotFile = pwobj.String(fnFasta)
                newItem._unitprotSeqLength = pwobj.Integer(sequenceLength(fnFasta))

            outputDatabaseID.append(newItem)

        fnAll = self._getPath("sequences.fasta")
        with open(fnAll, 'w') as outfile:
            for fname in fnList:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
                    outfile.write('\n\n')
        seqFile=ProteinSequenceFile()
        seqFile.setFileName(fnAll)

        self._defineOutputs(outputUniprot=outputDatabaseID)
        self._defineSourceRelation(self.inputListID, outputDatabaseID)
        self._defineOutputs(outputSequence=seqFile)
        self._defineSourceRelation(self.inputListID, seqFile)

    def _validate(self):
        errors=[]
        if not hasattr(self.inputListID.get().getFirstItem(),"_uniprotId"):
            errors.append("The set does not have an _uniprotId")
        return errors