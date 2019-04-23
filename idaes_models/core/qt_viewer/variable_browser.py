"""
A simple GUI viewer for Pyomo models.
"""
from __future__ import division
from __future__ import print_function
__author__ = "John Eslick"
__version__ = "1.0.0"

import warnings
try:
    from PyQt5 import QtCore
    from PyQt5.QtWidgets import QDialog, QAbstractItemView, QShortcut, QApplication
    from PyQt5 import QtGui
    from PyQt5 import uic
except:
    try:
        from PyQt4 import QtCore
        from PyQt4.QtGui import QDialog, QAbstractItemView, QShortcut, QApplication
        from PyQt4 import QtGui
        from PyQt4 import uic
    except:
        warnings.warn("Cannot import PyQt")


import os
from pyomo.environ import *
from pyomo.dae import *
import datetime


class VariableViewer(QDialog):
    def __init__(self, model=None):
        QDialog.__init__(self)
        path = os.path.dirname(__file__)
        self.ui = uic.loadUi(os.path.join(path, "variable_browser.ui"))
        self.model = model # Pyomo model
        print("Starting to Make Data Model")
        print(datetime.datetime.time(datetime.datetime.now()))
        datmodel = VariableDataModel(self, model)
        self.datmodel = datmodel
        self.ui.treeView.setModel(datmodel)
        #connect buttons
        self.ui.solve_button.setEnabled(False)
        self.ui.save_button.setEnabled(False)
        self.ui.load_button.setEnabled(False)
        self.ui.treeView.setSelectionBehavior(
            QAbstractItemView.SelectRows)
        self.ui.treeView.setSelectionMode(
            QAbstractItemView.ExtendedSelection)
        self.cp_shortcut = QShortcut(
            QtGui.QKeySequence("Ctrl+C"), self.ui)
        self.cp_shortcut.activated.connect(self.copyToClipboard)
        self.calc_shortcut = QShortcut(
            QtGui.QKeySequence("Ctrl+T"), self.ui)
        self.calc_shortcut.activated.connect(self.calculate)
        self.ui.show()

    def copyToClipboard(self):
        rows = set()
        cols = set()
        for i in self.ui.treeView.selectedIndexes():
            rows.add(i.row())
            cols.add(i.column())
        nrows = len(rows)
        ncols = len(cols)
        rows = list(sorted(rows))
        cols = list(sorted(cols))
        #make array to copy
        d = [None]*nrows
        for i in range(nrows): d[i] = [None]*ncols
        #get data
        for i in self.ui.treeView.selectedIndexes():
            e = i.data()
            if e is None:
                e = ""
            else:
                e = str(e)
            d[rows.index(i.row())][cols.index(i.column())] = e
        l = [None]*nrows
        for i in range(nrows):
            l[i] = "\t".join(d[i])
        s = "\n".join(l)
        clipboard = QApplication.clipboard()
        clipboard.setText(s)

    def calculate(self):
        for i in self.ui.treeView.selectedIndexes():
            if i.column() == 1:
                i.internalPointer().value(update=True, con=True)
            elif i.column() == 2:
                i.internalPointer().lb(update=True, con=True)
            elif i.column() == 3:
                i.internalPointer().ub(update=True, con=True)
        self.ui.update()

class TreeItem(object):
    def __init__(self, parent, row, o, otype=None):
        """
        Create a tree item
        """
        self.parent = parent
        self.row = row
        self.children = []
        self.data = o
        self.data_name = None
        self.data_value = None
        self.data_lb = None
        self.data_ub = None
        self.data_fixed = None
        self.data_active = None
        self.data_stale = None
        self.data_doc = None
        try:
            ut = o.unit_type[0]
        except:
            ut = None

        if ut is not None:
            self.data_type = ut
            try:
                self.data_active = o.active
            except:
                pass
        elif isinstance(o, pyomo.core.base.var.IndexedVar):
            self.data_type = 'IndexedVar'
        elif isinstance(o, DerivativeVar):
            self.data_type = 'DerivativeVar'
        elif isinstance(o, pyomo.core.base.var._VarData):
            self.data_type = 'Var'
            self.data_fixed = o.fixed
            self.data_stale = o.stale
        elif isinstance(o, pyomo.core.base.constraint.IndexedConstraint):
            self.data_type = 'IndexedConstraint'
        elif isinstance(o, pyomo.core.base.constraint._ConstraintData):
            self.data_type = 'Constraint'
        elif isinstance(o, Block) and not isinstance(o, pyomo.core.base.block._BlockData):
            self.data_type = 'IndexedBlock'
        elif isinstance(o, pyomo.core.base.block._BlockData):
            self.data_type = 'Block'
        else:
            self.data_type = str(type(o))
            try:
                self.data_fixed = o.fixed
            except:
                pass
        try:
            self.data_active = o.active
        except:
            pass
        try:
            self.data_doc = o.doc
        except:
            pass

    def name(self):
        if self.data_name is None:
            self.data_name = self.data.name
        return self.data_name

    def value(self, update=False, con=False):
        if self.data_value is not None and not update:
            return self.data_value
        if self.data_type == 'Var':
            try:
                self.data_value = value(self.data)
            except:
                pass
        elif self.data_type == 'Constraint' and con:
            try:
                self.data_value = value(self.data.body)
            except:
                pass
        return self.data_value

    def lb(self, update=False, con=False):
        if self.data_lb is not None and not update:
            return self.data_lb
        if self.data_type == 'Var':
            try:
                self.data_lb = self.data.lb
            except:
                pass
        elif self.data_type == 'Constraint' and con:
            try:
                self.data_lb = value(self.data.lower)
            except:
                pass
        return self.data_lb

    def ub(self, update=False, con=False):
        if self.data_ub is not None and not update:
            return self.data_ub
        if self.data_type == 'Var':
            try:
                self.data_ub = self.data.ub
            except:
                pass
        elif self.data_type == 'Constraint' and con:
            try:
                self.data_ub = value(self.data.upper)
            except:
                pass
        return self.data_ub

    def add_child(self, o):
        """
        Add a child tree item to this item
        """
        ti = TreeItem(self, 0, o)
        self.children.append(ti)
        return ti


class _ColIdx(object):
    NAME = 0
    VALUE = 1
    LB = 2
    UB = 3
    TYPE = 4
    FIXED = 5
    ACTIVE = 6
    STALE = 7
    DOC = 8
    N_CHILD = 9
    N_COL = 10


class VariableDataModel(QtCore.QAbstractItemModel):
    """
    This is a data model to provide the tree structure and information
    to the tree viewer
    """
    def __init__(self, parent, model):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self.col_head = [""]*_ColIdx.N_COL
        self.col_head[_ColIdx.NAME] = "Name"
        self.col_head[_ColIdx.VALUE] = "Value"
        self.col_head[_ColIdx.LB] = "L.B."
        self.col_head[_ColIdx.UB] = "U.B."
        self.col_head[_ColIdx.TYPE] = "Type"
        self.col_head[_ColIdx.FIXED] = "Fixed"
        self.col_head[_ColIdx.ACTIVE] = "Active"
        self.col_head[_ColIdx.STALE] = "Stale"
        self.col_head[_ColIdx.DOC] = "Doc"
        self.col_head[_ColIdx.N_CHILD] = "# Children"

        self.ncol = _ColIdx.N_COL
        self.model = model # Pyomo model
        self.rootItems = []
        self._create_tree()

    def _create_tree(self):
        """
        This contains a recursive function to go through a Pyomo model
        and create a tree structure with relevant information
        """
        def add_var(o, parent):
            """
            Add either a Var or an indexed/derivative Var element to
            the tree.
            """
            self._add_item(parent=parent, o=o)

        def add_constraint(o, parent):
            self._add_item(parent=parent, o=o)

        def o_to_tree(o, parent=None, count=None):
            """
            This is the recursive tree making function
            """
            if count is not None:
                count[0] += 1
                print(count[0])
            try:
                doc = o.doc
            except:
                doc = None
            if isinstance(o, Block) and not isinstance(o, pyomo.core.base.block._BlockData):
                p = self._add_item(parent=parent, o=o)
                for i, b in o.iteritems():
                    o_to_tree(b, p, count=count)
            elif isinstance(o, ConcreteModel):
                #treat the model the same at the block
                p = self._add_item(parent=parent, o=o)
                for no in o.component_objects(descend_into=False):
                    o_to_tree(no, p, count=count)
            elif isinstance(o, pyomo.core.base.block._BlockData):
                p = self._add_item(parent=parent, o=o)
                for no in o.component_objects(descend_into=False):
                    o_to_tree(no, p, count=count)
            elif isinstance(o, pyomo.core.base.var.IndexedVar):
                p = self._add_item(parent=parent, o=o)
                for i, no in o.iteritems():
                    if count is not None:
                        count[0] += 1
                        print(count[0])
                    add_var(no, p)
            elif isinstance(o, DerivativeVar):
                p = self._add_item(parent=parent, o=o)
                for i, no in o.iteritems():
                    if count is not None:
                        count[0] += 1
                        print(count[0])
                    add_var(no, p)
            elif isinstance(o, Var):
                add_var(o, parent)
            elif isinstance(o, pyomo.core.base.constraint.IndexedConstraint):
                p = self._add_item(parent=parent, o=o)
                for i, no in o.iteritems():
                    if count is not None:
                        count[0] += 1
                        print(count[0])
                    add_constraint(no, p)
            elif isinstance(o, Constraint):
                add_constraint(o, parent)

        #Enter the recursive function
        o_to_tree(self.model)
        print("Made Tree")
        print(datetime.datetime.time(datetime.datetime.now()))

    def _add_item(self, parent, o):
        """
        Add a root item in parent is None, otherwise add a child
        """
        if parent is None:
            ti = self._add_root_item(o)
        else:
            ti = parent.add_child(o)
        return ti

    def _add_root_item(self, o):
        """
        Add a root tree item
        """
        ti = TreeItem(None, 0, o)
        self.rootItems.append(ti)
        return ti

    def columnCount(self, parent=QtCore.QModelIndex()):
        """
        Return the number of columns
        """
        return self.ncol

    def rowCount(self, parent=QtCore.QModelIndex()):
        if not parent.isValid():
            return len(self.rootItems)
        return len(parent.internalPointer().children)

    def _o_value(self, o):
        if o:
            return "Not Calculated"
        elif self._is_var(o):
            return o.value
        else:
            return None

    def _o_doc(self, o):
        try:
            return o.doc
        except:
            return None

    def flags(self, index=QtCore.QModelIndex()):
        if index.column() in \
            [_ColIdx.VALUE,_ColIdx.LB,_ColIdx.UB,_ColIdx.FIXED, _ColIdx.ACTIVE]\
            and index.internalPointer().data_type == 'Var':
            return(QtCore.Qt.ItemIsEnabled |
                QtCore.Qt.ItemIsSelectable |
                QtCore.Qt.ItemIsEditable)
        elif index.column() == _ColIdx.ACTIVE and \
            index.internalPointer().data_active is not None:
            return(QtCore.Qt.ItemIsEnabled |
                QtCore.Qt.ItemIsSelectable |
                QtCore.Qt.ItemIsEditable)
        return(QtCore.Qt.ItemIsEnabled |
            QtCore.Qt.ItemIsSelectable)

    def setData(self, index, value, role=QtCore.Qt.EditRole):
        if role==QtCore.Qt.EditRole:
            if index.column() == _ColIdx.VALUE and \
                index.internalPointer().data_type == 'Var':
                try:
                    index.internalPointer().data.value = float(value)
                    index.internalPointer().data_value = float(value)
                except:
                    pass
            elif index.column() == _ColIdx.LB and \
                index.internalPointer().data_type == 'Var':
                try:
                    index.internalPointer().data.setlb(float(value))
                    index.internalPointer().data_lb = float(value)
                except:
                    pass
            elif index.column() == _ColIdx.UB and \
                index.internalPointer().data_type == 'Var':
                try:
                    index.internalPointer().data.setub(float(value))
                    index.internalPointer().data_ub = float(value)
                except:
                    pass
            elif index.column() == _ColIdx.FIXED and \
                index.internalPointer().data_type == 'Var':
                if value == 'true' or value == 'True' or value == "1":
                    index.internalPointer().data.fix()
                    index.internalPointer().data_fixed = True
                elif value == 'false' or value == 'False' or value == "0":
                    index.internalPointer().data.unfix()
                    index.internalPointer().data_fixed = False
            elif index.column() == _ColIdx.ACTIVE and \
                hasattr(index.internalPointer().data, "active"):
                if value == 'true' or value == 'True' or value == "1":
                    index.internalPointer().data.activate()
                    index.internalPointer().data_active = True
                elif value == 'false' or value == 'False' or value == "0":
                    index.internalPointer().data.deactivate()
                    index.internalPointer().data_active = False
        return True

    def data(self, index=QtCore.QModelIndex(), role=QtCore.Qt.DisplayRole):
        if role==QtCore.Qt.DisplayRole:
            if index.column() == _ColIdx.NAME:
                return index.internalPointer().name()
            elif index.column() == _ColIdx.VALUE:
                return index.internalPointer().value()
            elif index.column() == _ColIdx.LB:
                return index.internalPointer().lb()
            elif index.column() == _ColIdx.UB:
                return index.internalPointer().ub()
            elif index.column() == _ColIdx.TYPE:
                return index.internalPointer().data_type
            elif index.column() == _ColIdx.FIXED:
                return index.internalPointer().data_fixed
            elif index.column() == _ColIdx.ACTIVE:
                return index.internalPointer().data_active
            elif index.column() == _ColIdx.STALE:
                return index.internalPointer().data_stale
            elif index.column() == _ColIdx.DOC:
                return index.internalPointer().data_doc
            elif index.column() == _ColIdx.N_CHILD:
                n = len(index.internalPointer().children)
                if n == 0:
                    return None
                else:
                    return len(index.internalPointer().children)
            else:
                return None
        else:
            return

    def parent(self, index):
        if not index.isValid():
            return QtCore.QModelIndex()
        item = index.internalPointer()
        if item.parent is None:
            return QtCore.QModelIndex()
        else:
            return self.createIndex(0, 0, item.parent)

    def index(self, row, column, parent):
        if not parent.isValid():
            return self.createIndex(row, column, self.rootItems[row])
        parentItem = parent.internalPointer()
        return self.createIndex(row, column, parentItem.children[row])

    def headerData(self, i, orientation, role=QtCore.Qt.DisplayRole):
        """
        Return the column headings for the horizontal header and
        index numbers for the vertical header.
        """
        if orientation == QtCore.Qt.Horizontal and \
            role == QtCore.Qt.DisplayRole:
            return self.col_head[i]
        return None
