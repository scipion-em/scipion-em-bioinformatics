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

from math import ceil
import numpy as np

from pwem.protocols import EMProtocol
from pyworkflow.object import Float, Integer
from pyworkflow.protocol.params import PointerParam, EnumParam, StringParam, IntParam, FloatParam

class ProtBioinformaticsExportCSV(EMProtocol):
    """Export a set as a csv. It is located in the Run directory"""

    """
        Export metadata from a Scipion set into CSV format for external
        analysis and interoperability. The protocol converts the attributes
        associated with each object in the input set into a structured
        tabular file that can be opened in spreadsheet software or processed
        with data analysis tools.

        AI Generated:

        Export CSV (ProtBioinformaticsExportCSV) — User Manual
            Overview

            The Export CSV protocol converts the metadata stored in a Scipion
            EM set into a CSV text file. Its main purpose is to provide a
            simple and portable representation of the information associated
            with particles, volumes, micrographs, or other EM-related objects.

            In practical workflows, exporting metadata to CSV format is useful
            for downstream analysis, reporting, quality control, or integration
            with external software environments such as spreadsheets, statistical
            tools, or custom bioinformatics pipelines. The protocol focuses on
            exporting object attributes in a human-readable and machine-readable
            format.

            Inputs and General Workflow

            The protocol requires a single input set containing the objects
            whose metadata will be exported. During execution, the protocol
            iterates through all entries in the set and extracts their available
            attributes. The first processed object defines the column headers
            that will appear in the CSV file.

            Each object is written as a row in the output file, while every
            attribute is stored as a separate column. The resulting table
            preserves the structure of the metadata in a format that can be
            easily reused outside Scipion.

            Outputs and Interpretation

            After execution, the protocol generates an output CSV file inside
            the protocol run directory. The file contains the exported metadata
            organized in tabular form, with one row per object and one column
            per attribute.

            This output can be used for data inspection, statistical analysis,
            dataset filtering, or interoperability with third-party software.
            Since the protocol exports raw attribute values directly, the
            interpretation of the data depends on the type of objects contained
            in the original set.

            Practical Recommendations

            The protocol is especially useful when users need to exchange
            metadata with external collaborators, generate reports, or perform
            custom analyses outside the Scipion ecosystem. Before exporting,
            it is advisable to verify that the input set contains consistent
            and meaningful metadata fields, since these attributes will define
            the structure of the resulting CSV table.

            Final Perspective

            The Export CSV protocol provides a lightweight and practical way
            to transform Scipion metadata into a universally accessible format.
            By simplifying interoperability and external analysis, it helps
            integrate cryo-EM workflows with broader computational and data
            science environments.
        """
    _label = 'export csv'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam, pointerClass="EMSet",
                       label='Set:', allowsNull=False)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('exportStep')

    def exportStep(self):
        fh = open(self._getPath("output.csv"),"w")
        lineNo = 0
        lineHdr = ""
        for entry in self.inputSet.get():
            line = ""
            for key, value in entry.getAttributes():
                if lineNo == 0:
                    if lineHdr!="":
                        lineHdr+="; "
                    lineHdr += key
                if line!="":
                    line += "; "
                line+=str(value.get())
            if lineNo==0:
                fh.write(lineHdr+"\n")
            fh.write(line+"\n")
            lineNo+=1
        fh.close()
