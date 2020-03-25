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

from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import AtomicStructHandler
from pyworkflow.protocol.params import (EnumParam, PointerParam, StringParam)

class ProtAtomStructDali(EMProtocol):
    """Query Dali server (http://ekhidna2.biocenter.helsinki.fi/dali) with a structure.
    The server returns other structures that are structurally similar"""
    methodsDict = {0: 'search', 1: 'pdb25'}
    _label = 'dali'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, pointerClass="AtomStruct",
                       label='Atomic Structure:', allowsNull=False,
                       help="Input atomic structure for the query")
        form.addParam('method', EnumParam,
                         choices=[val for key, val in sorted(self.methodsDict.items())],
                         label="Method:", default=0,
                         help="Select query mode (see http://ekhidna2.biocenter.helsinki.fi/dali)")
        form.addParam('title', StringParam,
                         label="Job title:", default="Scipion search",
                         help="Put a label that helps you to identify the job")
        form.addParam('email', StringParam,
                         label="Email:", default="",
                         help="The web page will send an email to this address when the calculations are finished. "
                              "Copy the URL at the email in the Analyze Results of this protocol")

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('searchStep',self.inputStructure.get().getFileName())

    def searchStep(self, structFileName):
        outFileName = self._getExtraPath("atomStruct.pdb")
        aStruct1 = AtomicStructHandler(structFileName)
        aStruct1.write(outFileName)
        args='-F "file1=@%s" -F "method=%s" -F "title=%s"  -F "address=%s" http://ekhidna.biocenter.helsinki.fi/cgi-bin/dali/dump.cgi' %\
             (outFileName,self.methodsDict[self.method.get()],self.title.get(),self.email.get())
        self.runJob("curl",args)

    # --------------------------- UTILS functions ------------------
    def _validate(self):
        errors = []
        if self.email.get()=="":
            errors.append("The email cannot be empty")
        from pyworkflow.utils.which import which
        if which("curl") == "":
            errors.append("Cannot find curl in the path. Install with apt-get install curl, yum install curl or equivalent in your system")
        if which("wget") == "":
            errors.append("Cannot find curl in the path. Install with apt-get install wget, yum install wget or equivalent in your system")
        return errors

    def _summary(self):
        summary = []
        return summary

    def _citations(self):
        return ['Holm2019']
