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
    debugger.HandleCommand(
        'type summary add --summary-string "${var.m_data.__elems_[0]}" *real<.*>'
    )
    # lupnt Vector and Matrix
    debugger.HandleCommand(
        "type synthetic add -x lupnt::(Vector|Matrix)* --python-class eigenlldb.MatrixChildProvider"
    )
    debugger.HandleCommand(
        "type summary add -x lupnt::(Vector|Matrix)* --python-function eigenlldb.summary"
    )


def summary(valobj, internal_dict, options):
    names = [str(x.GetName()) for x in valobj.get_value_child_list()]
    values = [float(x.GetValue()) for x in valobj.get_value_child_list()]

    if len(names[0]) == 3:
        # Vector "[i]"
        row_idxs = [int(x[1:-1]) for x in names]
        rows = max(row_idxs) + 1
        col_idxs = [0] * rows
        cols = 1
    else:
        # Matrix "[i,j]"
        tmp = [x[1:-1].split(",") for x in names]
        row_idxs = [int(x[0]) for x in tmp]
        rows = max(row_idxs) + 1
        col_idxs = [int(x[1]) for x in tmp]
        cols = max(col_idxs) + 1

    # Array
    array = [[0.0] * cols for _ in range(rows)]
    for i, j, x in zip(row_idxs, col_idxs, values):
        array[i][j] = x

    return print_aligned(array)


def print_aligned(matrix):
    num_columns = len(matrix[0])
    max_before_decimal = [0] * num_columns
    max_after_decimal = [0] * num_columns

    # Step 1: Find the maximum length before and after decimal for each column
    for row in matrix:
        for i, value in enumerate(row):
            parts = str(value).split(".")
            max_before_decimal[i] = max(max_before_decimal[i], len(parts[0]))
            if len(parts) > 1:
                max_after_decimal[i] = max(max_after_decimal[i], len(parts[1]))

    # Step 2: Format each number with padding and accumulate in a string
    txt = "\n["
    for r, row in enumerate(matrix):
        line = "[" if r == 0 else " ["
        for i, value in enumerate(row):
            value = float(format_element(value))
            before, after = (
                str(value).split(".") if "." in str(value) else (str(value), "")
            )
            padding_left = " " * (max_before_decimal[i] - len(before))
            padding_right = " " * (max_after_decimal[i] - len(after))
            formatted_value = f"{padding_left}{before}.{after}{padding_right}"
            line += formatted_value + (" " if i < num_columns - 1 else "")
        txt += line + "]\n"

    return txt


def format_element(x):
    if x == 0.0:
        return "0.0"
    else:
        return f"{x}"


def format_array(array):
    txt = "\n"
    txt += "["
    for i, row in enumerate(array):
        if i != 0:
            txt += "\n "
        txt += "["
        txt += " ".join(format_element(x) for x in row)
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
