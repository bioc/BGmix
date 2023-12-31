#  This file is part of BGX, the Bayesian Gene eXpression program.
#  Copyright 2007 Ernest Turro <ernest.turro@ic.ac.uk>
#
#  BGX is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License, version 2, as
#  published by the Free Software Foundation.
#
#  BGX is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


# Initialise.
AC_INIT([BGmix], [1.1.1], [a.m.lewin@ic.ac.uk])
AC_CONFIG_SRCDIR([src/BGmix_main.cpp])
#AC_CONFIG_HEADERS([src/config.h])

# Checks for programs.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CXX=`"${R_HOME}/bin/R" CMD config CXX`


BGMIX_CFLAGS="$CFLAGS"
BGMIX_CXXFLAGS="$CXXFLAGS"

# We are compiling for R.
BGMIX_CXXFLAGS="$BGMIX_CXXFLAGS -DUSING_R"

# Get C99 flag(s) (if any), and append to BGX_CFLAGS.
PRIORCC=$CC
AC_PROG_CC_C99 # May change CC
for i in $CC; do 
  if test `expr "$PRIORCC" : ".*$i"` -eq 0;then 
    BGMIX_CFLAGS="$BGMIX_CFLAGS $i"
  fi;
done

#Link against math library (sets LIBS)
AC_CHECK_LIB([m],[cos],,[AC_MSG_ERROR([Please install math library])])

# Checks for Boost libraries.
#AX_BOOST([1.31.0], BGMIX_CPPFLAGS="$BOOST_CPPFLAGS")
if test \( -f ./boostIncl.tar.gz \) -a \( ! -d ./boostIncl \);
  then
    echo "Untarring Boost header files..."
    tar -xzf  ./boostIncl.tar.gz
fi;
BGMIX_CPPFLAGS="-I ../boostIncl"

# Checks for typedefs, structures, and compiler characteristics.
#AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE

# Finish.
AC_SUBST(BGMIX_CPPFLAGS)
AC_SUBST(BGMIX_CXXFLAGS)
AC_SUBST(BGMIX_CFLAGS)
AC_SUBST(LIBS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

