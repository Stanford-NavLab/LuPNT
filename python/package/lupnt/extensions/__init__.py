# import lupnt.pylupnt as lupnt


# # List of all lupnt C++ number types exported to Python. Extend this as other
# # number types are exported!
# _lupnt_number_types = [
#     lupnt.real1st,
#     lupnt.real2nd,
#     lupnt.real3rd,
#     lupnt.real4th,
# ]


# # Define the __format__ methods for all lupnt number types. This is needed so
# # that we can write formatted strings such as f"{x:.3f}" and avoid a Python
# # runtime error.
# def _lupnt_number_format(self, spec):
#     return format(self.val(), spec)


# for numbertype in _lupnt_number_types:
#     numbertype.__format__ = _lupnt_number_format

# # -------------------------------------------------------------------------------------
# # The code below is needed so that lupnt numbers can be used in plotly figures.
# # It teaches PlotlyJSONEncoder on how to encode an lupnt number to JSON.
# try:
#     from plotly.utils import PlotlyJSONEncoder

#     plotly_default = PlotlyJSONEncoder.default

#     def is_lupnt_number_type(o):
#         return any(t for t in _lupnt_number_types if isinstance(o, t))

#     def plotly_default_new(self, o):
#         if is_lupnt_number_type(o):
#             return o.val()
#         else:
#             return plotly_default(self, o)

#     PlotlyJSONEncoder.default = plotly_default_new
# except:
#     pass
