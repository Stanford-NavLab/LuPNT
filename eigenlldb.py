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

import numpy as np

MAX_ROWS = 10
MAX_COLS = 10


def __lldb_init_module(debugger, internal_dict):
    # Real
    debugger.HandleCommand(
        'type summary add --summary-string "${var.m_data.__elems_[0]}" lupnt::Real'
    )
    debugger.HandleCommand(
        'type summary add --summary-string "${var.m_data.__elems_[0]}" autodiff::detail::Real<1, double>'
    )
    # lupnt Vec and Mat
    debugger.HandleCommand(
        "type synthetic add -x lupnt::(Vec|Mat)* --python-class eigenlldb.MatChildProvider"
    )
    debugger.HandleCommand(
        "type summary add -x lupnt::(Vec|Mat)* --python-function eigenlldb.summary"
    )
    debugger.HandleCommand(
        "type synthetic add -x Eigen::(Vector|Matrix)* --python-class eigenlldb.MatChildProvider"
    )
    debugger.HandleCommand(
        "type summary add -x Eigen::(Vector|Matrix)* --python-function eigenlldb.summary"
    )


def summary(valobj, internalr_dict, options):
    names = [str(x.GetName()) for x in valobj.get_value_child_list()]
    values = [float(x.GetValue()) for x in valobj.get_value_child_list()]
    matrix = MatChildProvider(valobj, internalr_dict)
    filetered_names = [x for x in names if x not in ("m_rows", "m_cols")]
    filtered_values = [
        v for n, v in zip(names, values) if n not in ("m_rows", "m_cols")
    ]

    cols = (
        matrix._cols_compile_time
        if matrix._cols_compile_time != -1
        else int(values[names.index("m_cols")])
    )
    rows = (
        matrix._rows_compile_time
        if matrix._rows_compile_time != -1
        else int(values[names.index("m_rows")])
    )
    cols_desc = "dynamic" if matrix._cols_compile_time == -1 else "static"
    rows_desc = "dynamic" if matrix._rows_compile_time == -1 else "static"

    if len(names[-1]) == 3:
        # Vec "[i]"
        row_idxs = [int(x[1:-1]) for x in filetered_names]
        n_rows = max(row_idxs) + 1
        col_idxs = [0] * n_rows
        n_cols = 1
    else:
        # Mat "[i,j]"
        tmp = [x[1:-1].split(",") for x in filetered_names]
        row_idxs = [int(x[0]) for x in tmp]
        n_rows = max(row_idxs) + 1
        col_idxs = [int(x[1]) for x in tmp]
        n_cols = max(col_idxs) + 1

    # Array
    array = [[0.0] * n_cols for _ in range(n_rows)]
    for i, j, x in zip(row_idxs, col_idxs, filtered_values):
        array[i][j] = x

    summary_txt = "(" + str(rows) + "," + str(cols) + ") "
    summary_txt += "(" + cols_desc + "," + rows_desc + ")\n"
    summary_txt += print_aligned(
        array, all_rows=rows == n_rows, all_columns=cols == n_cols
    )
    # summary_txt += "\n" + str(np.array(array))
    return summary_txt


def print_aligned(matrix, all_rows=True, all_columns=True):
    num_columns = len(matrix[0])
    max_before_decimal = [0] * num_columns
    max_after_decimal = [0] * num_columns

    # Step 1: Find the maximum length before and after decimal for each column
    for row in matrix:
        for i, value in enumerate(row):
            parts = format_element(value).split(".")
            max_before_decimal[i] = max(max_before_decimal[i], len(parts[0]))
            if len(parts) > 1:
                max_after_decimal[i] = max(max_after_decimal[i], len(parts[1]))

    # Step 2: Format each number with padding and accumulate in a string
    txt = "["
    for r, row in enumerate(matrix):
        line = "[" if r == 0 else " ["
        for i, value in enumerate(row):
            value = format_element(value)
            before, after = (
                str(value).split(".") if "." in str(value) else (str(value), "")
            )
            padding_left = " " * (max_before_decimal[i] - len(before))
            padding_right = " " * (max_after_decimal[i] - len(after))
            formatted_value = f"{padding_left}{before}.{after}{padding_right}"
            line += formatted_value + (" " if i < num_columns - 1 else "")
        if not all_columns:
            line += " ...]"
        else:
            line += "]"

        txt += line
        if r < len(matrix) - 1:
            txt += "\n"
    if not all_rows:
        txt += "\n [" + " " * (len(line) - 6) + "...]"
    txt += "]"
    return txt


def format_element(x):
    if x == 0.0:
        return "0.0"
    else:
        return f"{x:.6g}"


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


class MatChildProvider:
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

        self._fixed_rows = max_rows != -1
        self._fixed_cols = max_cols != -1
        self._fixed_storage = self._fixed_rows and self._fixed_cols

    def num_children(self):
        return (
            min(self._cols(), MAX_COLS) * min(self._rows(), MAX_ROWS)
            + (not self._fixed_rows)
            + (not self._fixed_cols)
        )

    def get_child_index(self, name):
        return int(name.lstrip("[").rstrip("]"))

    def get_child_at_index(self, index):
        storage = self._valobj.GetChildMemberWithName("m_storage")
        data = storage.GetChildMemberWithName("m_data")

        if not self._fixed_rows:
            if index == 0:
                return storage.GetChildMemberWithName("m_rows")
            index -= 1
        if not self._fixed_cols:
            if (index == 0 and not self._fixed_rows) or (
                index == 1 and self._fixed_rows
            ):
                return storage.GetChildMemberWithName("m_cols")
            index -= 1
        if self._fixed_storage:
            data = data.GetChildMemberWithName("array")

        n_cols = min(self._cols(), MAX_COLS)
        n_rows = min(self._rows(), MAX_ROWS)
        if self._row_major:
            row = index // n_cols
            col = index % n_cols
            new_index = row * self._cols() + col
        else:
            row = index % n_rows
            col = index // n_rows
            new_index = col * self._rows() + row
        offset = self._scalar_size * new_index

        if n_cols == 1:
            name = "[{}]".format(row)
        elif n_rows == 1:
            name = "[{}]".format(col)
        else:
            name = "[{},{}]".format(row, col)

        if str(self._scalar_type) in ("double", "int", "bool"):
            # double/int/bool
            child = data.CreateChildAtOffset(name, offset, self._scalar_type)
        else:
            # Real
            if self._fixed_storage:
                child = data.GetChildAtIndex(new_index)
            else:
                child = data.CreateChildAtOffset(name, offset, self._scalar_type)
            child = child.GetChildMemberWithName("m_data")
            child = child.GetChildMemberWithName("__elems_")
            offset = child.GetChildAtIndex(0).GetByteSize() * new_index
            scalar_type = child.GetChildAtIndex(0).GetType()
            child = child.CreateChildAtOffset(name, 0, scalar_type)
        return child

    def _cols(self):
        if self._cols_compile_time == -1:
            storage = self._valobj.GetChildMemberWithName("m_storage")
            cols = storage.GetChildMemberWithName("m_cols")
            n_cols = cols.GetValueAsUnsigned()
        else:
            n_cols = self._cols_compile_time
        return n_cols

    def _rows(self):
        if self._rows_compile_time == -1:
            storage = self._valobj.GetChildMemberWithName("m_storage")
            rows = storage.GetChildMemberWithName("m_rows")
            n_rows = rows.GetValueAsUnsigned()
        else:
            n_rows = self._rows_compile_time
        return n_rows
