# -*- coding: utf-8 -*-
# This file is part of Eigen, a lightweight C++ template library
# for linear algebra.
#
# Copyright (C) 2021 Huang, Zhaoquan <zhaoquan2008@hotmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Pretty printers for Eigen::Matrix to use with LLDB debugger
#
# Usage:
# 1. Add the following line (change it according to the path to this file)
#    to the file ~/.lldbinit (create one if it doesn't exist):
#        `command script import /path/to/eigenlldb.py`
# 2. Inspect the variables in LLDB command line
#        `frame variable`

import lldb
from typing import List
import bisect


def __lldb_init_module(debugger, internal_dict):
    # real
    debugger.HandleCommand(
        'type summary add --summary-string "${var.m_data._M_elems[0]}" real'  # ???
    )
    # lupnt Vector and Matrix
    debugger.HandleCommand(
        'type summary add --summary-string "${var.m_data.__elems_[0]}" *real<.*>'
    )
    debugger.HandleCommand(
        "type summary add -x lupnt::(Vector|Matrix)* --python-function eigenlldb.summary"
    )
    debugger.HandleCommand(
        "type synthetic add -x lupnt::(Vector|Matrix)* --python-class eigenlldb.MatrixChildProvider"
    )


def summary(valobj, internal_dict, options):
    """
    Create summary for Eigen
    """
    valtype = valobj.GetType().GetCanonicalType()
    scalar_type = valtype.GetTemplateArgumentType(0)
    if not scalar_type.IsValid():
        # In the case that scalar_type is invalid on LLDB 9.0 on Windows with CLion
        storage = valobj.GetChildMemberWithName("m_storage")
        data = storage.GetChildMemberWithName("m_data")
        data_type = data.GetType()
        if data_type.IsPointerType():
            scalar_type = data.GetType().GetPointeeType()
        else:
            scalar_type = (
                data.GetChildMemberWithName("array").GetType().GetArrayElementType()
            )
    scalar_size = scalar_type.GetByteSize()
    name = valtype.GetName()
    # template_begin = name.find("<")
    # template_end = name.find(">")
    # template_args = name[(template_begin + 1) : template_end].split(",")

    # find the last '>' to support nested template
    # split template arguments without splitting nested template
    # example
    # in  : "Eigen::Matrix<autodiff::detail::Real<1, double>, -1, 1, 0, -1, 1>""
    # out : ["autodiff::detail::Real<1, double>", "-1", "1", "0", "-1", "1"]
    template_args = extract_args(name)

    rows_compile_time = int(template_args[1])
    cols_compile_time = int(template_args[2])
    row_major = (int(template_args[3]) & 1) != 0

    max_rows = int(template_args[4])
    max_cols = int(template_args[5])
    fixed_storage = max_rows != -1 and max_cols != -1
    storage = valobj.GetChildMemberWithName("m_storage")
    data = storage.GetChildMemberWithName("m_data")

    if cols_compile_time == -1:
        # Dynamic size
        storage = valobj.GetChildMemberWithName("m_storage")
        cols = storage.GetChildMemberWithName("m_cols")
        cols = cols.GetValueAsUnsigned()
    else:
        # Fixed size
        cols = cols_compile_time

    if rows_compile_time == -1:
        # Dynamic size
        storage = valobj.GetChildMemberWithName("m_storage")
        rows = storage.GetChildMemberWithName("m_rows")
        rows = rows.GetValueAsUnsigned()
    else:
        # Fixed size
        rows = rows_compile_time

    if fixed_storage:
        data_ = data.GetChildMemberWithName("array")
    else:
        data_ = data

    array = [[] for _ in range(rows)]
    for index in range(rows * cols):
        offset = scalar_size * index

        if row_major:
            row = index // cols
            col = index % cols
        else:
            row = index % rows
            col = index // rows

        if cols == 1:
            name = "[{}]".format(row)
        elif rows == 1:
            name = "[{}]".format(col)
        else:
            name = "[{},{}]".format(row, col)

        if str(scalar_type) == "double":
            # double
            child = data_.CreateChildAtOffset(name, offset, scalar_type)
            val = float(child.GetValue())
        else:
            # real
            child = data_.CreateChildAtOffset(name, offset, scalar_type)
            child = child.GetChildMemberWithName("m_data")
            child = child.GetChildMemberWithName("__elems_")  # __elems_/_M_elems
            val = float(child.GetChildAtIndex(0).GetValue())
            der = float(child.GetChildAtIndex(1).GetValue())

        array[row].append(val)

    return format_array(array)


def format_elem(x):
    if x == 0.0:
        return "0.0"
    else:
        return f"{x:.4f}"


def format_array(array):
    txt = "\n"
    txt += "["
    for i, row in enumerate(array):
        if i != 0:
            txt += "\n "
        txt += "["
        txt += " ".join(format_elem(x) for x in row)
        txt += "]"
    txt += "]"
    return txt


def split_template_args(s):
    result = []
    brace_count = 0
    current_arg = ""

    for char in s:
        if char == "<":
            brace_count += 1
        elif char == ">":
            brace_count -= 1
        elif char == "," and brace_count == 0:
            result.append(current_arg.strip())
            current_arg = ""
            continue

        current_arg += char

    if current_arg:
        result.append(current_arg.strip())

    return result


def extract_args(s):
    s = s.strip().split("<", 1)[-1].rsplit(">", 1)[0]
    return split_template_args(s)


class MatrixChildProvider:
    _valobj: lldb.SBValue
    _scalar_type: lldb.SBType
    _scalar_size: int
    _rows_compile_time: int
    _cols_compile_time: int
    _row_major: bool
    _fixed_storage: bool

    def __init__(self, valobj, internal_dict):
        self._valobj = valobj
        valtype = valobj.GetType().GetCanonicalType()

        scalar_type = valtype.GetTemplateArgumentType(0)
        if not scalar_type.IsValid():
            # In the case that scalar_type is invalid on LLDB 9.0 on Windows with CLion
            storage = valobj.GetChildMemberWithName("m_storage")
            data = storage.GetChildMemberWithName("m_data")
            data_type = data.GetType()
            if data_type.IsPointerType():
                scalar_type = data.GetType().GetPointeeType()
            else:
                scalar_type = (
                    data.GetChildMemberWithName("array").GetType().GetArrayElementType()
                )
        self._scalar_type = scalar_type
        self._scalar_size = self._scalar_type.GetByteSize()

        name = valtype.GetName()
        template_args = extract_args(name)

        self._rows_compile_time = int(template_args[1])
        self._cols_compile_time = int(template_args[2])
        self._row_major = (int(template_args[3]) & 1) != 0

        max_rows = int(template_args[4])
        max_cols = int(template_args[5])
        self._fixed_storage = max_rows != -1 and max_cols != -1

    def num_children(self):
        return self._cols() * self._rows()

    def get_child_index(self, name):
        pass

    def get_child_at_index(self, index):
        storage = self._valobj.GetChildMemberWithName("m_storage")
        data = storage.GetChildMemberWithName("m_data")
        offset = self._scalar_size * index

        if self._row_major:
            row = index // self._cols()
            col = index % self._cols()
        else:
            row = index % self._rows()
            col = index // self._rows()
        if self._fixed_storage:
            data = data.GetChildMemberWithName("array")
        if self._cols() == 1:
            name = "[{}]".format(row)
        elif self._rows() == 1:
            name = "[{}]".format(col)
        else:
            name = "[{},{}]".format(row, col)

        if str(self._scalar_type) == "double":
            # double
            child = data.CreateChildAtOffset(name, offset, self._scalar_type)
        else:
            # real
            child = data.CreateChildAtOffset(name, offset, self._scalar_type)
            child = child.GetChildMemberWithName("m_data")
            child = child.GetChildMemberWithName("__elems_")
            scalar_type = child.GetChildAtIndex(0).GetType()
            child = child.CreateChildAtOffset(name, offset, scalar_type)
        return child

    def _cols(self):
        if self._cols_compile_time == -1:
            storage = self._valobj.GetChildMemberWithName("m_storage")
            cols = storage.GetChildMemberWithName("m_cols")
            return cols.GetValueAsUnsigned()
        else:
            return self._cols_compile_time

    def _rows(self):
        if self._rows_compile_time == -1:
            storage = self._valobj.GetChildMemberWithName("m_storage")
            rows = storage.GetChildMemberWithName("m_rows")
            return rows.GetValueAsUnsigned()
        else:
            return self._rows_compile_time
