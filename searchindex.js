Search.setIndex({"docnames": ["PDEs/PDEs", "PDEs/PDEs-Poisson", "PDEs/PDEs-flows", "PDEs/PDEs-iterative", "PDEs/PDEs-waves", "README", "basiclinalg/bla", "basiclinalg/bla-docu", "basiclinalg/bla-expressions", "basiclinalg/bla-git", "basiclinalg/bla-libs", "basiclinalg/bla-performance", "basiclinalg/bla-python", "basiclinalg/bla-refreshCpp", "intro", "markdown", "markdown-notebooks", "notebooks"], "filenames": ["PDEs/PDEs.md", "PDEs/PDEs-Poisson.md", "PDEs/PDEs-flows.md", "PDEs/PDEs-iterative.md", "PDEs/PDEs-waves.md", "README.md", "basiclinalg/bla.md", "basiclinalg/bla-docu.md", "basiclinalg/bla-expressions.md", "basiclinalg/bla-git.md", "basiclinalg/bla-libs.md", "basiclinalg/bla-performance.md", "basiclinalg/bla-python.md", "basiclinalg/bla-refreshCpp.md", "intro.md", "markdown.md", "markdown-notebooks.md", "notebooks.ipynb"], "titles": ["<span class=\"section-number\">9. </span>Partial differential equations", "<span class=\"section-number\">10. </span>Poisson Equation and the 5-point stencil", "<span class=\"section-number\">13. </span>Navier-Stokes Equation", "<span class=\"section-number\">11. </span>Iterative Solvers", "<span class=\"section-number\">12. </span>Wave equations", "IntroSC", "<span class=\"section-number\">1. </span>Overview", "<span class=\"section-number\">4. </span>Creating documentation", "<span class=\"section-number\">5. </span>Expression templates", "<span class=\"section-number\">2. </span>Clone a project from github", "<span class=\"section-number\">7. </span>Overview of libraries", "<span class=\"section-number\">8. </span>Performance", "<span class=\"section-number\">6. </span>Python bindings", "<span class=\"section-number\">3. </span>Refresh your C++", "Introduction to Scientific Computing", "Markdown Files", "Notebooks with MyST Markdown", "Content with notebooks"], "terms": {"poisson": [0, 14], "5": [0, 14, 17], "point": [0, 14], "stencil": [0, 14], "iter": [0, 14], "solver": [0, 14], "wave": [0, 14], "navier": [0, 14], "stoke": [0, 14], "lectur": [5, 14], "note": [5, 15], "introduct": 5, "scientif": 5, "comput": [5, 8, 9], "read": 5, "onlin": [5, 8], "we": [6, 7, 8, 9, 13, 14], "ar": [6, 8, 13, 15, 16], "develop": [6, 9, 14], "librari": [6, 8, 9, 14], "oper": [6, 8, 9, 13], "For": [6, 13, 15, 17], "thi": [6, 8, 13, 14, 15, 16, 17], "implement": [6, 8], "class": [6, 8, 9, 13, 14], "vector": [6, 8, 9, 13, 14], "matric": [6, 14], "modern": 6, "c": [6, 8, 9, 14], "make": [6, 7, 9, 17], "them": 6, "avail": [6, 9], "python": [6, 14], "simpl": [6, 15], "usabl": 6, "practic": 6, "softwar": [6, 9, 14], "team": 6, "us": [6, 7, 8, 9, 13, 15, 16], "git": [6, 8, 9], "start": [6, 15, 16], "maintain": 6, "document": [6, 14, 15, 16, 17], "from": [6, 8, 13, 14, 17], "begin": [6, 17], "learn": [6, 14], "trick": [6, 8], "like": [6, 8, 9, 15, 16], "express": [6, 14], "templat": [6, 13, 14], "combin": [6, 8], "perform": [6, 8, 14], "readabl": 6, "work": [8, 17], "refresh": 14, "your": [7, 9, 14, 15, 16, 17], "creat": [8, 13, 14, 17], "bind": 14, "overview": [14, 15], "readthedoc": 7, "base": [7, 8, 9, 16], "sphinx": [7, 15], "quick": 7, "step": [7, 9], "setup": 7, "mkdir": [7, 9], "doc": 7, "cd": [7, 9], "quickstart": 7, "conf": 7, "py": 7, "you": [7, 8, 9, 13, 15, 16, 17], "can": [7, 8, 9, 13, 15, 16, 17], "select": 7, "theme": 7, "exampl": [7, 8, 13, 15, 17], "html_theme": 7, "sphinx_rtd_them": 7, "regist": 7, "github": [7, 14], "account": 7, "On": 7, "dash": 7, "board": 7, "see": [7, 8, 13, 15, 16, 17], "public": [7, 8, 13], "project": [7, 13, 14], "import": [7, 8, 17], "think": 8, "code": [8, 9, 13, 15, 16], "z": 8, "x": 8, "3": 8, "y": 8, "where": 8, "what": [8, 13], "happen": [8, 13], "first": 8, "doubl": [8, 13], "call": [8, 13, 15], "gener": 8, "temporari": [8, 13], "new": [8, 13], "sum": [8, 13], "final": 8, "later": 8, "copi": 8, "memori": 8, "write": [8, 15, 16], "an": [8, 13, 14, 15], "old": 8, "style": [8, 13], "function": [8, 9, 13, 15], "size_t": [8, 13], "i": [8, 13, 14], "0": [8, 13, 17], "size": [8, 13], "seem": 8, "much": 8, "more": [8, 9, 16, 17], "effici": 8, "object": [8, 13], "result": 8, "valu": [8, 13], "store": [8, 15], "right": 8, "posit": 8, "immedi": 8, "provid": [8, 9], "wai": 8, "eleg": 8, "version": [8, 9], "explicit": 8, "loop": 8, "The": [8, 9, 13, 15, 16], "doe": 8, "return": [8, 13], "logic": 8, "inform": [8, 16, 17], "am": 8, "second": 8, "scale": 8, "scalar": 8, "Such": 8, "built": [8, 16], "compil": 8, "someth": 8, "type": [8, 13], "sumexpr": 8, "scaleexpr": 8, "encod": 8, "paradigm": 8, "One": [8, 9], "also": [8, 15, 16, 17], "refer": [8, 15], "lazi": 8, "evalu": 8, "sinc": 8, "just": [8, 15], "when": [8, 15, 16], "written": [8, 15, 16], "destin": 8, "asc": 8, "bla": 8, "clone": [8, 14], "branch": 8, "expr": 8, "here": [8, 13, 15, 17], "have": [8, 13, 14, 16], "all": [8, 13, 15, 16], "matrix": [8, 9, 13], "typenam": [8, 13], "t": [8, 13, 15, 17], "auto": 8, "view": 8, "const": [8, 13], "static_cast": 8, "ta": 8, "tb": 8, "a_": 8, "b_": 8, "b": [8, 13], "fanci": 8, "deriv": 8, "give": 8, "itself": 8, "argument": [8, 13], "so": [8, 16], "static": 8, "up": 8, "cast": 8, "idiom": 8, "known": 8, "curious": 8, "recur": 8, "pattern": 8, "crtp": 8, "barton": 8, "nackman": 8, "In": [8, 14, 15], "break": 8, "todd": 8, "veldhuizen": 8, "wa": 8, "introduc": [8, 13], "howev": 8, "back": [8, 9], "1995": 8, "too": 8, "technolog": 8, "If": [8, 9, 13, 16], "upcast": 8, "both": [8, 15], "its": [8, 13], "member": [8, 9, 13], "But": [8, 17], "realli": 8, "good": [8, 13], "nest": 8, "ye": 8, "It": [8, 9, 15], "inlin": [8, 15], "whole": 8, "flow": 8, "data": [8, 13, 17], "optim": 8, "everyth": 8, "singl": 8, "To": 8, "verifi": 8, "look": [8, 9], "assembl": 8, "There": [8, 17], "tool": [8, 15], "explor": 8, "sourc": [8, 9], "displai": [8, 16], "allow": [8, 15], "choos": 8, "between": 8, "lot": [8, 15, 17], "flag": 8, "popular": 9, "control": 9, "system": 9, "free": 9, "open": 9, "platform": 9, "cloud": 9, "servic": 9, "prerequisit": 9, "instal": 9, "environ": 9, "e": [9, 14], "g": 9, "visual": 9, "studio": 9, "within": [9, 13], "vs": 9, "cmake": 9, "draft": 9, "either": 9, "vscode": 9, "window": 9, "now": 9, "configur": 9, "build": [9, 15], "find": 9, "file": [9, 16], "cmakelist": 9, "txt": 9, "describ": 9, "structur": [9, 15], "again": 9, "do": [9, 14, 15, 17], "termin": 9, "go": 9, "directori": 9, "run": [9, 16], "should": [9, 16], "abl": 9, "execut": [9, 16], "add": 9, "want": [9, 17], "form": [9, 13], "2": [9, 16], "4": [9, 16, 17], "student": [9, 14], "fork": 9, "name": [9, 13], "invit": 9, "colleagu": 9, "repeat": [9, 13], "befor": [9, 13], "homework": 9, "follow": [9, 13, 15, 16], "product": 9, "transpos": 9, "some": [9, 13, 15, 17], "test": 9, "commit": 9, "push": 9, "chang": 9, "command": [9, 16], "line": [9, 15, 16, 17], "instruct": [9, 16], "cc": 9, "h": 9, "m": 9, "ad": [9, 13], "other": [9, 16], "fetch": 9, "merg": 9, "contribut": 9, "experi": 9, "same": [9, 13, 15], "region": 9, "assum": 13, "basic": [6, 13], "knowledg": 13, "onli": 13, "shortli": 13, "featur": 13, "program": 13, "languag": [13, 15], "list": 13, "book": [13, 15, 16, 17], "size_": 13, "data_": 13, "v": 13, "nullptr": 13, "std": 13, "swap": 13, "delet": 13, "v2": 13, "defin": [13, 16], "element": 13, "given": 13, "one": [13, 15], "precis": 13, "entri": 13, "A": 13, "constructor": 13, "destroi": 13, "destructor": 13, "move": 13, "semant": 13, "11": 13, "initi": 13, "which": [13, 14, 16], "anywai": 13, "steal": 13, "resourc": 13, "rule": 13, "five": 13, "privat": 13, "mean": [13, 17], "invis": 13, "outsid": 13, "interact": [13, 17], "via": 13, "method": 13, "assign": 13, "chosen": [13, 14], "instead": 13, "bracket": 13, "consist": 13, "two": [13, 14, 15, 16], "recommend": 13, "coher": 13, "includ": [13, 16, 17], "convent": 13, "variabl": 13, "indent": 13, "space": 13, "after": 13, "try": 13, "googl": 13, "guid": [13, 17], "joachim": 14, "sch\u00f6berl": 14, "tu": 14, "wien": 14, "institut": 14, "analysi": 14, "mani": [14, 15, 16], "challeng": 14, "requir": 14, "non": 14, "trivial": 14, "solut": 14, "tackel": 14, "research": 14, "divers": 14, "skill": 14, "approach": 14, "typic": 14, "applic": 14, "area": 14, "numer": 14, "partial": 14, "differenti": 14, "equat": 14, "design": [6, 14], "appli": 14, "mathemat": 14, "earli": 14, "stage": 14, "whether": 15, "s": [15, 16, 17], "content": [15, 16], "jupyt": [15, 16, 17], "notebook": 15, "ipynb": 15, "regular": 15, "md": [15, 16], "ll": 15, "flavor": 15, "help": 15, "get": [15, 16], "show": [15, 16], "off": [15, 16], "syntax": 15, "stand": 15, "markedli": 15, "text": [15, 16], "slight": 15, "variat": 15, "commonmark": 15, "small": 15, "extens": 15, "ecosystem": 15, "about": [15, 16, 17], "most": 15, "power": 15, "thei": 15, "kind": 15, "markup": 15, "serv": 15, "similar": 15, "purpos": 15, "wherea": 15, "span": 15, "accept": 15, "differ": 15, "input": 15, "those": 15, "depend": 15, "specif": 15, "being": 15, "render": 15, "special": 15, "box": 15, "cite": 15, "bibtex": 15, "holdgraf_evidence_2014": 15, "hdhpk14": 15, "moreov": 15, "insert": 15, "bibliographi": 15, "page": [15, 16], "must": 15, "properli": 15, "bib": 15, "christoph": 15, "ramsai": 15, "holdgraf": 15, "wendi": 15, "de": 15, "heer": 15, "brian": 15, "n": [15, 17], "paslei": 15, "robert": 15, "knight": 15, "evid": 15, "predict": 15, "human": 15, "auditori": 15, "cortex": 15, "intern": 15, "confer": 15, "cognit": 15, "neurosci": 15, "brisban": 15, "australia": 15, "2014": 15, "frontier": 15, "starter": 15, "jupyterbook": 15, "org": 15, "let": 16, "detail": 16, "With": 16, "direct": 16, "print": 16, "ani": 16, "block": 16, "default": 16, "kernel": 16, "output": 16, "rest": 16, "jupytext": 16, "convert": 16, "support": 16, "thing": 16, "need": 16, "understand": 16, "how": [6, 16], "top": 16, "presenc": 16, "That": 16, "d": 16, "treat": 16, "init": 16, "path": 16, "markdownfil": 16, "As": 17, "emb": 17, "imag": 17, "html": 17, "etc": 17, "post": 17, "add_": 17, "math": 17, "align": 17, "mbox": 17, "la_": 17, "tex": 17, "end": 17, "sure": 17, "escap": 17, "dollar": 17, "sign": 17, "keep": 17, "well": 17, "check": 17, "out": 17, "sampl": 17, "matplotlib": 17, "rcparam": 17, "cycler": 17, "pyplot": 17, "plt": 17, "numpi": 17, "np": 17, "ion": 17, "contextlib": 17, "exitstack": 17, "0x120050b50": 17, "fix": 17, "random": 17, "state": 17, "reproduc": 17, "seed": 17, "19680801": 17, "10": 17, "logspac": 17, "1": 17, "100": 17, "randn": 17, "ii": 17, "rang": 17, "arrai": 17, "cmap": 17, "cm": 17, "coolwarm": 17, "ax": 17, "prop_cycl": 17, "color": 17, "linspac": 17, "line2d": 17, "custom_lin": 17, "lw": 17, "fig": 17, "subplot": 17, "figsiz": 17, "plot": 17, "legend": 17, "cold": 17, "medium": 17, "hot": 17, "linear": 6, "algebra": 6, "cpu": 6, "algorithm": 6, "util": 6, "potenti": 6}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"partial": 0, "differenti": 0, "equat": [0, 1, 2, 4], "poisson": 1, "5": 1, "point": 1, "stencil": 1, "navier": 2, "stoke": 2, "iter": 3, "solver": 3, "wave": 4, "introsc": 5, "basic": 14, "linear": 14, "algebra": 14, "creat": [7, 16], "document": 7, "express": 8, "templat": 8, "clone": 9, "project": 9, "from": 9, "github": 9, "extend": 9, "asc": 9, "bla": 9, "team": 9, "work": 9, "overview": [6, 10], "librari": 10, "perform": 11, "python": 12, "bind": 12, "refresh": 13, "your": 13, "c": 13, "introduct": 14, "scientif": 14, "comput": 14, "pde": 14, "markdown": [15, 16, 17], "file": 15, "what": 15, "myst": [15, 16, 17], "sampl": 15, "role": 15, "direct": 15, "citat": 15, "learn": 15, "more": 15, "notebook": [16, 17], "an": 16, "exampl": 16, "cell": 16, "quickli": 16, "add": 16, "yaml": 16, "metadata": 16, "content": 17, "code": 17, "block": 17, "output": 17}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinxcontrib.bibtex": 9, "sphinx": 56}})