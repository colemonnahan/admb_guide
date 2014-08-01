(TeX-add-style-hook "admb_guide"
 (lambda ()
    (LaTeX-add-bibliographies)
    (LaTeX-add-labels
     "sec:outfiles"
     "sec:restart"
     "sec:diag"
     "sec:startvals"
     "sec:MH"
     "tab:mh_args"
     "fig:mcrb"
     "sec:hybrid"
     "eq:motion")
    (TeX-run-style-hooks
     "float"
     "graphicx"
     "amssymb"
     "amsmath"
     "geometry"
     "margin=1in"
     "latex2e"
     "art10"
     "article"
     "")))

