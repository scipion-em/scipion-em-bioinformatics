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

from pyworkflow.protocol.params import (StringParam, PointerParam)
from pwem.protocols import EMProtocol
from bioinformatics.utils.utils import *

class ProtBioinformaticsRaptorX(EMProtocol):
    """This is a wrapper to http://raptorx.uchicago.edu/StructPredV2/predict/. The underlying program
       predicts the 3D structure of a protein from its aminoacid sequence"""

    """
        Submit protein amino acid sequences to the RaptorX web server for
        automated protein structure prediction. The protocol acts as a
        wrapper around the RaptorX Structure Prediction service, allowing
        users to send FASTA sequences directly from Scipion workflows and
        integrate remote structural prediction into bioinformatics and
        structural biology analyses.

        AI Generated:

        RaptorX Structure Prediction (ProtBioinformaticsRaptorX) — User Manual
            Overview

            The RaptorX Structure Prediction protocol provides an interface
            to the RaptorX online prediction server for estimating the
            three-dimensional structure of proteins from amino acid
            sequences. Its primary purpose is to simplify the submission
            of protein sequences to the external prediction platform while
            maintaining integration with Scipion bioinformatics workflows.

            In structural biology and computational modeling workflows,
            this protocol is especially useful when experimental structures
            are unavailable and predicted models are required for further
            interpretation, docking studies, comparative analyses, or
            exploratory structural investigations.

            Inputs and General Workflow

            The protocol requires a protein sequence provided as a FASTA
            entry or sequence object. During execution, the sequence is
            extracted and submitted to the RaptorX prediction server using
            a command-line HTTP request. Users may optionally define a job
            title and an email address to facilitate identification of the
            prediction task and receive notification once the remote
            calculations are completed.

            The protocol automatically generates a default execution label
            when no custom title is provided, ensuring traceability of
            prediction jobs within workflow environments.

            Remote Prediction Submission

            Internally, the protocol relies on external command-line tools
            such as curl to communicate with the RaptorX web service. The
            sequence is transmitted through an HTTP request directed to
            the prediction endpoint hosted by the RaptorX infrastructure.

            Since the prediction process is executed remotely, the protocol
            itself does not compute structural models locally. Instead, it
            serves as a submission and integration layer connecting Scipion
            workflows with the external prediction platform.

            Validation and Execution Requirements

            Before execution, the protocol validates that the input
            contains a proper FASTA sequence and verifies the availability
            of required command-line utilities such as curl and wget in
            the local system environment.

            Stable internet connectivity is essential because the protocol
            depends entirely on communication with the remote RaptorX
            server infrastructure.

            Final Perspective

            The RaptorX Structure Prediction protocol provides a simple and
            integrated mechanism for connecting protein sequence analysis
            workflows with remote structure prediction services. By
            automating sequence submission and workflow integration, the
            protocol facilitates structural modeling studies when
            experimentally determined protein structures are not available.
        """
    _label = 'raptorX'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSeq', PointerParam, pointerClass="SetOfDatabaseID, Sequence",
                       label='Sequence:', allowsNull=False,
                       help="Input FASTA sequence, the set must have a single structure")
        form.addParam('title', StringParam,
                         label="Job title:", default="",
                         help="Put a label that helps you to identify the job")
        form.addParam('email', StringParam,
                         label="Email:", default="",
                         help="The web page will send an email to this address when the calculations are finished. "
                              "Copy the URL at the email in the Analyze Results of this protocol")

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('searchStep')

    def searchStep(self):
        sequence = copyFastaSequenceAndRead(self)

        title = self.title.get()
        if not title:
            title="ScipionRun%s"%os.path.split(self._getPath())[1][0:6]

        args='-F "jobname=%s" -F "email=%s" -F "seqeunces=%s" http://raptorx.uchicago.edu/StructPredV2/predict/' %\
             (title, self.email.get(), sequence)
        self.runJob("curl",args)


    # --------------------------- UTILS functions ------------------
    def _validate(self):
        errors = checkInputHasFasta(self)
        from pyworkflow.utils.which import which
        if which("curl") == "":
            errors.append("Cannot find curl in the path. Install with apt-get install curl, yum install wget or equivalent in your system")
        if which("wget") == "":
            errors.append("Cannot find curl in the path. Install with apt-get install wget, yum install wget or equivalent in your system")
        return errors

    def _citations(self):
        return ['Zheng2019'] # ***
