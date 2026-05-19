# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *              Alberto M. Parra Pérez (amparraperez@gmail.com)
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

import glob
import os


from pwem.protocols import EMProtocol
import pyworkflow.object as pwobj
from pyworkflow.utils.path import copyFile
from pyworkflow.protocol.params import PathParam, StringParam, BooleanParam
from bioinformatics.objects import SmallMolecule, SetOfSmallMolecules
from bioinformatics import Plugin

class ProtBioinformaticsImportSmallMolecules(EMProtocol):
    """Import small molecules from a directory. Each molecule should be in a separate file.
       Smiles (.smi), SDF (.sdf), Tripos Mol2 (.mol2), Maestro (.mae, .maegz), PDB blocks (.pdb)
       are accepted.
    """
    """
        Import and organize small-molecule structures from different chemical
        file formats into a Scipion-compatible dataset. The protocol supports
        individual molecular files as well as multi-entry chemical databases,
        automatically converting and separating compounds into independent
        molecular objects for downstream analysis.

        AI Generated:

        Import Small Molecules (ProtBioinformaticsImportSmallMolecules) — User Manual
            Overview

            The Import Small Molecules protocol provides a flexible framework
            for incorporating chemical compound libraries into Scipion workflows.
            Its main purpose is to convert collections of molecular structures
            into a standardized internal representation that can later be used
            for visualization, docking, screening, or bioinformatics analyses.

            The protocol supports several commonly used molecular formats,
            including SMILES, SDF, Mol2, Maestro, and PDB files. These formats
            are widely employed in cheminformatics, structural biology, and
            virtual screening workflows, allowing the protocol to integrate
            heterogeneous molecular datasets into a unified environment.

            Inputs and General Workflow

            The protocol allows users to import molecules either from a directory
            containing independent molecular files or from a single multi-entry
            file containing multiple compounds. Input selection can be controlled
            through file patterns and wildcard expressions, facilitating the
            import of large chemical libraries distributed across different
            folders or naming conventions.

            When a multi-entry file is provided, the protocol automatically
            separates individual molecules into independent files. Multi-molecule
            Mol2 and SDF datasets are parsed entry by entry, while CSV-based
            SMILES collections are converted into separate molecular structure
            files identified by their compound identifiers.

            After processing the input data, the protocol generates a Scipion
            set of small molecules in which every imported compound is stored
            as an individual molecular object ready for downstream processing.

            Molecular Representation and Visualization

            During import, the protocol optionally generates molecular preview
            images using RDKit utilities. These images provide a rapid visual
            representation of the imported ligands and facilitate inspection
            of chemical datasets directly within the workflow environment.

            Image generation is automatically limited for very large datasets
            in order to avoid excessive computational overhead. Unsupported
            or failed visualizations are safely handled without interrupting
            the import process.

            Outputs and Interpretation

            After execution, the protocol produces a SetOfSmallMolecules object
            containing all successfully imported compounds. Each molecule retains
            its associated structural file and, when possible, a generated
            molecular image preview.

            The resulting dataset can be directly reused in subsequent workflows
            involving docking, ligand analysis, structural comparison, or
            bioinformatics pipelines. Since the protocol preserves the original
            molecular information, the biological interpretation depends on the
            quality and consistency of the imported compound library.

            Practical Recommendations

            The protocol is particularly useful when preparing ligand libraries
            for virtual screening or structural analysis workflows. For optimal
            organization, it is recommended to use standardized identifiers and
            consistent molecular formats before import.

            Multi-entry SDF and Mol2 databases should be carefully validated
            prior to execution to ensure that compound delimiters and identifiers
            are correctly formatted. When working with extremely large datasets,
            users should also consider the computational cost associated with
            molecular image generation.

            Final Perspective

            The Import Small Molecules protocol provides a practical bridge
            between external chemical databases and Scipion molecular workflows.
            By supporting multiple chemical formats and automated dataset
            organization, the protocol simplifies the integration of small-
            molecule information into structural biology and computational
            bioinformatics environments.
        """
    _label = 'import small mols'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('multiple', BooleanParam, default=True,
                      label='Each file is a molecule')
        form.addParam('filesPath', PathParam, condition='multiple',
                      label="Files directory",
                      help="Directory with the files you want to import.\n\n"
                           "The path can also contain wildcards to select"
                           "from several folders. \n\n"
                           "Examples:\n"
                           "  ~/project/data/day??_files/\n"
                           "Each '?' represents one unknown character\n\n"
                           "  ~/project/data/day*_files/\n"
                           "'*' represents any number of unknown characters\n\n"
                           "  ~/project/data/day##_files/\n"
                           "'##' represents two digits that will be used as "
                           "file ID\n\n"
                           "NOTE: wildcard characters ('*', '?', '#') "
                           "cannot appear in the actual path.)")
        form.addParam('filesPattern', StringParam,  condition='multiple',
                      label='Pattern',
                      default="*",
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.\n\n"
                           "NOTE: wildcards and special characters "
                           "('*', '?', '#', ':', '%') cannot appear in the "
                           "actual path.\n\n"
                           "You may create small molecules from Smiles (.smi), Tripos Mol2 (.mol2), "
                           "SDF (.sdf), Maestro (.mae, .maegz), or PDB blocks (.pdb)")

        form.addParam('filePath', PathParam, condition='not multiple',
                      label='File', help='Allowed formats: \n '
                                         ' - CSV smiles (ID, compound; this is downloaded from ZINC) \n'
                                         ' - Mol2 (Multiple mol2 file) \n'
                                         ' - SDF (Multiple sdf file)')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importStep')

    def importStep(self):

        if not self.multiple.get():

            fnSmall = self._getExtraPath(os.path.split(self.filePath.get())[1])
            copyFile(self.filePath.get(), fnSmall)

            if fnSmall.endswith(".mol2"): # Multiple mol2
                with open(fnSmall) as f:
                    lines = f.readlines()
                    i = 0
                    for line in lines:
                        if line.startswith("#"):
                            continue

                        if line.startswith("@<TRIPOS>MOLECULE"):
                            zincID = lines[i+1].split()[0]
                            fname_small = self._getExtraPath("%s.mol2" %zincID)
                            f_small = open(fname_small, 'w+')
                            f_small.write(line)

                        f_small.write(line)
                        try:
                            if lines[i+1].startswith("@<TRIPOS>MOLECULE"):
                                f_small.close()
                        except:
                            f_small.close()

                        i += 1

                os.remove(fnSmall)


            elif fnSmall.endswith(".sdf"): # Multiple sdf
                with open(fnSmall) as f:
                    lines = f.readlines()
                    lines2write = []
                    i = 0
                    for line in lines:
                        lines2write.append(i)

                        if line.startswith("$$$$"):
                            zincID = lines[i-5].split()[0]
                            fname_small = self._getExtraPath("%s.sdf" % zincID)
                            f_small = open(fname_small, 'w+')
                            for l in lines2write:
                                f_small.write(lines[l])

                            lines2write = []
                            f_small.close()

                        i += 1

                os.remove(fnSmall)

            else:
                fh = open(self.filePath.get())
                for line in fh.readlines():
                    tokens = line.split(',')
                    if len(tokens)==2:
                        fhSmile = open(self._getExtraPath(tokens[0].strip()+".smi"),'w')
                        fhSmile.write(tokens[1].strip()+"\n")
                        fhSmile.close()
                fh.close()

        else:
            for filename in glob.glob(os.path.join(self.filesPath.get(), self.filesPattern.get())):
                fnSmall = self._getExtraPath(os.path.split(filename)[1])
                copyFile(filename, fnSmall)

        outputSmallMolecules = SetOfSmallMolecules().create(path=self._getPath(),suffix='SmallMols')


        for fnSmall in glob.glob(self._getExtraPath("*")):
            smallMolecule = SmallMolecule(smallMolFilename=fnSmall)

            if len(os.listdir(self._getExtraPath())) <= 100: # costly
                if not fnSmall.endswith('.mae') or not fnSmall.endswith('.maegz'):
                    fnRoot = os.path.splitext(os.path.split(fnSmall)[1])[0]
                    fnOut = self._getExtraPath("%s.png" % fnRoot)
                    args = Plugin.getPluginHome('utils/rdkitUtils.py') + " draw %s %s" % (fnSmall, fnOut)
                    try:
                        Plugin.runRDKit(self, "python3", args)
                        smallMolecule._PDBLigandImage = pwobj.String(fnOut)
                    except:
                        smallMolecule._PDBLigandImage = pwobj.String("Not available")

            outputSmallMolecules.append(smallMolecule)
        self._defineOutputs(outputSmallMols=outputSmallMolecules)
