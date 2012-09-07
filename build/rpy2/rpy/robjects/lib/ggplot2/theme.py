from rpy2.robjects.packages import importr

ggplot2_pack = importr('ggplot2')


class Options(robjects.Vector):
   def __init__(self, obj):
      self.__sexp__ = obj.__sexp__

   def __repr__(self):
      s = '<instance of %s : %i>' %(type(self), id(self)) 
      return s


class Theme(Options):
   pass

class Blank(Theme):
    _constructor = ggplot2_pack.theme_blank
    @classmethod
    def new(cls):
        res = cls(cls._constructor())
        return res


class Grey(Theme):
    _constructor = ggplot2_pack.theme_grey
    @classmethod
    def new(cls, base_size = 12):
       res = cls(cls._constructor(base_size = base_size))
       return res


class Rect(Theme):
    _constructor = ggplot2_pack.theme_rect
    @classmethod
    def new(cls, fill = robjects.NA_Logical, colour = "black", 
            size = 0.5, linetype = 1):
       res = cls(cls._constructor(fill = fill, colour = colour, 
                                  size = size, linetype = linetype))
       return res


class Segment(Theme):
    _constructor = ggplot2_pack.theme_rect
    @classmethod
    def new(cls, colour = 'black', size = 0.5, linetype = 1):
       res = cls(cls._constructor(colour = colour, size = size,
                                  linetype = linetype))
       return res


# Theme text is not a vector :/
class Text(robjects.Function):
    _constructor = ggplot2_pack.theme_text
    @classmethod
    def new(cls, family = "", face = "plain", colour = "black", size = 10,
            hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1.1):
       res = cls(cls._constructor(family = family, face = face, 
                                  colour = colour, size = size,
                                  hjust = hjust, vjust = vjust, 
                                  angle = angle, lineheight = lineheight))
       return res


class BW(Theme):
    _constructor = ggplot2_pack.theme_bw
    @classmethod
    def new(cls, base_size = 12):
       res = cls(cls._constructor(base_size = base_size))
       return res


class Gray(Theme):
    _constructor = ggplot2_pack.theme_gray
    @classmethod
    def new(cls, base_size = 12):
       res = cls(cls._constructor(base_size = base_size))
       return res


class Line(Theme):
    _constructor = ggplot2_pack.theme_line
    @classmethod
    def new(cls, colour = 'black', size = 0.5, linetype = 1):
       res = cls(cls._constructor(colour = colour, size = size,
                                  linetype = linetype))
       return res
