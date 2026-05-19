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

class ProtBioinformaticsListOperate(EMProtocol):
    """Filter a set by a column value or keep just a few columns"""

    """
        Perform filtering, sorting, counting, intersection, and column
        selection operations on Scipion metadata sets. The protocol
        provides a flexible framework for manipulating EM-related
        datasets according to user-defined criteria and attribute-based
        conditions.

        AI Generated:

        Operate Set (ProtBioinformaticsListOperate) — User Manual
            Overview

            The Operate Set protocol allows users to manipulate metadata
            collections stored in Scipion EM sets through a variety of
            operations such as filtering, sorting, selecting columns,
            identifying unique entries, counting occurrences, or finding
            intersections between datasets. Its main purpose is to provide
            a practical and flexible mechanism for reorganizing and refining
            metadata during cryo-EM and bioinformatics workflows.

            In practical data analysis pipelines, metadata operations are
            essential for selecting biologically relevant subsets,
            simplifying large datasets, preparing data for downstream
            analysis, or generating curated collections for visualization
            and interpretation. The protocol focuses on attribute-based
            manipulation while preserving compatibility with the original
            Scipion object structure.

            Inputs and General Workflow

            The protocol requires an input EM set containing the metadata
            entries to be processed. Depending on the selected operation,
            users may additionally provide filtering conditions, numerical
            thresholds, sorting directions, lists of columns to preserve,
            or a second reference set for intersection analysis.

            During execution, the protocol iterates through the entries
            contained in the input set and applies the selected operation
            directly to the object attributes. New output entries are then
            generated while preserving the original object type and
            metadata structure.

            Filtering and Attribute-Based Selection

            One of the central functionalities of the protocol is metadata
            filtering based on attribute values. Users can define numerical
            or textual conditions including equality, inequality, interval
            comparisons, substring matching, and prefix or suffix searches.
            This enables efficient extraction of biologically meaningful
            subsets such as particles above a confidence threshold,
            structures containing specific annotations, or entries matching
            experimental conditions.

            From a practical perspective, filtering operations are useful
            for quality control, dataset refinement, and preparation of
            curated subsets prior to structural interpretation or machine
            learning analysis. Careful selection of filter conditions is
            important to avoid excluding biologically relevant entries.

            Dataset Reduction and Structural Organization

            The protocol also supports operations aimed at simplifying or
            reorganizing metadata collections. Users may keep only selected
            columns, generate unique entries according to a chosen
            attribute, or sort datasets in ascending or descending order.
            Additional operations allow extraction of the highest or lowest
            scoring entries based on absolute counts or percentiles.

            These functionalities are particularly useful when handling
            large cryo-EM datasets where reducing redundancy, prioritizing
            top-ranked results, or reorganizing metadata improves workflow
            efficiency and interpretability.

            Counting and Intersection Analysis

            The protocol includes counting operations that determine how
            many times specific attribute values appear within a dataset.
            This allows users to identify repeated classifications,
            frequently occurring annotations, or dominant structural
            categories within large collections.

            Intersection analysis enables comparison between two metadata
            sets by preserving only entries whose identifiers are shared
            between datasets. This functionality is useful for cross-study
            validation, consensus analysis, or integration of independent
            computational workflows.

            Outputs and Interpretation

            After execution, the protocol generates a new EM set containing
            the processed metadata entries. The output preserves the
            original object structure while reflecting the selected
            transformations or filtering criteria.

            Depending on the operation performed, the resulting dataset may
            represent a refined subset, a reordered collection, a reduced
            metadata representation, or an annotated set containing
            additional counting information. The interpretation of the
            output therefore depends directly on the selected operation and
            the biological context of the original data.

            Practical Recommendations

            The protocol is especially useful when preparing datasets for
            downstream structural analysis, machine learning workflows,
            statistical evaluation, or publication-ready curation. Before
            applying restrictive filters or percentile-based selections, it
            is advisable to verify the consistency and biological meaning
            of the selected metadata attributes.

            In large-scale cryo-EM workflows, combining filtering,
            uniqueness detection, sorting, and intersection operations can
            substantially improve dataset organization and reproducibility.

            Final Perspective

            The Operate Set protocol provides a versatile framework for
            metadata manipulation and dataset refinement within the
            Scipion ecosystem. By enabling flexible attribute-based
            operations on EM sets, the protocol facilitates data curation,
            interoperability, and biologically meaningful analysis across
            structural biology workflows.
        """
    _label = 'operate set'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('operation', EnumParam, choices=['Filter', 'Keep columns', 'Unique', 'Top N', 'Bottom N',
                                                       'Top %', 'Bottom %', 'Count', 'Intersection', 'Sort'],
                      label='Operation', default=0,
                      help='In intersection, we keep those entries of the Set to filter whose identifier (filter column) '
                           'are in the second set')
        form.addParam('inputSet', PointerParam, pointerClass="EMSet",
                       label='Set to filter:', allowsNull=False)
        form.addParam('secondSet', PointerParam, pointerClass="EMSet", condition="operation==8", # Intersection between 2 Sets
                       label='Set with IDs:', allowsNull=True)
        form.addParam('filterColumn', StringParam,
                       label='Filter column:', condition='(operation!=1)',
                       help='It must exist in the input object.')
        form.addParam('filterOp', EnumParam, choices=['==', '>', '>=', '<', '<=', '!=', 'between', 'startswith',
                                                      'endswith', 'contains', 'does not startwith',
                                                      'does not end with', 'does not contain'],
                       label='Filter operation:', condition='(operation==0)')
        form.addParam('filterValue', StringParam,
                       label='Value:', condition='(operation==0)',
                       help='Value to use in the filter')
        form.addParam('filterValue2', StringParam,
                       label='Lower Value:', condition='(operation==0 and filterOp==6)',
                       help='Value to use in the filter')
        form.addParam('keepColumns', StringParam,
                       label='Keep columns:', condition='(operation==1)',
                       help='They must exist in the input database list. Separated by semicolons (e.g. column1 ; column2 ; ...)')
        form.addParam('N', IntParam,
                       label='N:', default=500, condition='(operation==3 or operation==4)')
        form.addParam('percentile', FloatParam,
                       label='Percentile:', default=5, condition='(operation==5 or operation==6)',
                       help='Between 0 and 100')
        form.addParam('direction', EnumParam, choices=['Ascending', 'Descending'], default=0,
                       label='Sorting direction:', condition='(operation==9)')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('operateStep')

    def operateStep(self):

        outputSet = self.inputSet.get().create(self._getPath())

        if self.operation.get()==0:
            # Filter columns
            referenceValue = self.filterValue.get()
            value = self.inputSet.get().getFirstItem().getAttributeValue(self.filterColumn.get())

            if isinstance(value,float):
                referenceValue = float(referenceValue)
            elif isinstance(value,int):
                referenceValue = int(referenceValue)

            filterOp = self.filterOp.get()
            if filterOp == 6:
                referenceValue2 = self.filterValue2.get()

                if isinstance(value, float):
                    referenceValue2 = float(referenceValue2)
                elif isinstance(value, int):
                    referenceValue2 = int(referenceValue2)

            for oldEntry in self.inputSet.get():
                value = oldEntry.getAttributeValue(self.filterColumn.get())
                if isinstance(value, Float):
                    value = float(value)
                elif isinstance(value, Integer):
                    value = int(value)


                add = False
                if filterOp == 0: # ==
                    add = value==referenceValue
                elif filterOp == 1: # >
                    add = value>referenceValue
                elif filterOp == 2:  # >=
                    add = value > referenceValue
                elif filterOp == 3:  # <
                    add = value < referenceValue
                elif filterOp == 4:  # <=
                    add = value <= referenceValue
                elif filterOp == 5:  # !=
                    add = value != referenceValue
                elif filterOp == 6:  # between
                    add = (value <= referenceValue and value >= referenceValue2)
                elif filterOp == 7:  #startswith
                    add = value.startswith(referenceValue)
                elif filterOp == 8:  # endswith
                    add = value.endswith(referenceValue)
                elif filterOp == 9:  # contains
                    add = referenceValue in value
                elif filterOp == 10:  # does not startswith
                    add = not (value.startswith(referenceValue))
                elif filterOp == 11:  # does not endswith
                    add = not (value.endswith(referenceValue))
                elif filterOp == 12:  # does not contains
                    add = not (referenceValue in value)
                if add:
                    newEntry = self.inputSet.get().ITEM_TYPE()
                    newEntry.copy(oldEntry)
                    outputSet.append(newEntry)

        elif self.operation.get()==1:
            # Keep columns
            keepList=[x.strip() for x in self.keepColumns.get().split()]

            ignoreList=[]
            for name, _ in self.inputSet.get().getFirstItem().getAttributes():
                if not name in keepList:
                    ignoreList.append(name)
            for oldEntry in self.inputSet.get():
                newEntry = self.inputSet.get().ITEM_TYPE()
                newEntry.copy(oldEntry,ignoreAttrs=ignoreList)
                outputSet.append(newEntry)

        elif self.operation.get()==2:
            # Unique
            found={}
            for oldEntry in self.inputSet.get():
                value = oldEntry.getAttributeValue(self.filterColumn.get())
                if not value in found:
                    found[value] = True
                    newEntry = self.inputSet.get().ITEM_TYPE()
                    newEntry.copy(oldEntry)
                    outputSet.append(newEntry)

        elif self.operation.get()>=3 and self.operation.get()<=6:
            # Top N, Bottom N,Top %, Bottom %
            V = []
            for entry in self.inputSet.get():
                V.append(entry.getAttributeValue(self.filterColumn.get()))
            V.sort()
            op = self.operation.get()
            if op==3:
                threshold = V[-self.N.get()]
            elif op==4:
                threshold = V[self.N.get()-1]
            elif op==5:
                threshold = V[-ceil(self.percentile.get()/100*len(V))]
            elif op==6:
                threshold = V[ceil(self.percentile.get()/100*len(V))-1]

            for oldEntry in self.inputSet.get():
                value = oldEntry.getAttributeValue(self.filterColumn.get())
                if (op==3 or op==5) and value>=threshold:
                    newEntry = self.inputSet.get().ITEM_TYPE()
                    newEntry.copy(oldEntry)
                    outputSet.append(newEntry)
                elif (op==4 or op==6) and value<=threshold:
                    newEntry = self.inputSet.get().ITEM_TYPE()
                    newEntry.copy(oldEntry)
                    outputSet.append(newEntry)

        elif self.operation.get()==7:
            # Count the number of entries that are the same
            count={}
            for oldEntry in self.inputSet.get():
                value = oldEntry.getAttributeValue(self.filterColumn.get())
                if not value in count:
                    count[value] = 0
                count[value] += 1

            for oldEntry in self.inputSet.get():
                value = oldEntry.getAttributeValue(self.filterColumn.get())

                newEntry = self.inputSet.get().ITEM_TYPE()
                newEntry.copy(oldEntry)
                newEntry.count = Integer(count[value])
                outputSet.append(newEntry)

        elif self.operation.get()==8:
            # Intersection between 2 Sets
            secondSet={}
            for entry in self.secondSet.get():
                value = entry.getAttributeValue(self.filterColumn.get())
                if not value in secondSet:
                    secondSet[value] = True

            for oldEntry in self.inputSet.get():
                value = oldEntry.getAttributeValue(self.filterColumn.get())

                if value in secondSet:
                    newEntry = self.inputSet.get().ITEM_TYPE()
                    newEntry.copy(oldEntry)
                    outputSet.append(newEntry)

        elif self.operation.get()==9:
            # Sort
            V = []
            newEntries = []
            for entry in self.inputSet.get():
                V.append(entry.getAttributeValue(self.filterColumn.get()))
                newEntry = self.inputSet.get().ITEM_TYPE()
                newEntry.copy(entry)
                newEntry.cleanObjId()
                newEntries.append(newEntry)
            if self.direction.get()==0:
                idxSort = np.argsort(V)
            else:
                idxSort = np.argsort(-np.asarray(V))

            for idx in idxSort:
                outputSet.append(newEntries[idx])

        if len(outputSet)>0:
            self._defineOutputs(output=outputSet)
            self._defineSourceRelation(self.inputSet, outputSet)
