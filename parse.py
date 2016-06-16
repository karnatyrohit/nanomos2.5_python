
import numpy as np

def parse(fin, p):
    #
    # function [p] = parser(fin)
    #
    # September 2000
    #
    # Author: Jeffery L. Gray, grayj@purdue.edu
    #         School of Electrical and Computer Engineering
    #         Purdue University, West Lafayette, IN 47907
    #__________________________
    #
    # Input -
    #
    #       fin: file id of input file containting statements
    #
    # Output -
    #
    #       p: a structure, defined below
    #
    # The following variables are always returned:
    #
    # p.ncard - is a string containing the name of the input card
    # p.err - error messages, as follows:
    #    normal return conditions:
    #         p.err=-1; p.errmess='end of file';
    #         p.err=1; p.errmess='no errors detected on assignment statement';
    #         p.err=2; p.errmess='no errors detected on title statement';
    #    error detected:
    #         p.err=999; p.errmess=strcat('cannot mix numbers and strings
    #in: ',p.var(i).name);
    #
    # parser.m reads stements from the file fin which are of the form:
    #
    #       *titlecard  anything at all
    #
    #   or
    #
    #       cardname  var1=value1,var2=string2 var3=value3
    #       +     array1=va1/va1/va3, array2=str1/str2/str3/str4
    #       +     var4 value2
    #
    # Lines beginning with a blank or a $ are ignored.
    #
    # Variables have 3 types: number, string, or empty
    #
    # Statements may be any number of lines. The card name
    # must start in column 1. The continuation symbol, +,
    # must also appear in column 1.
    #
    # Commas or blanks are assumed to be separaters. Any
    # number of separaters may appear between assignments.
    # An assignment cannot contain any blanks, i.e.
    #      block = - 12.0
    # is not valid. It must be written as
    #      block=-12.0   instead.
    #
    # Multiple values can be assigned to the variable
    # by separating the values by /'s. For example
    # wl=.3/.4/.5/.6/.7/.8
    #    and
    # ohmic=yes/no/no/no/yes
    #
    #
    # Examples-
    #
    #       |<-- column 1
    #
    #####   newcard  x=3.0,rat=cat  alpha=1.0e12
    #       +  wl=.3/.4/.5/inf/1e500 string=ab/cde/fghi/j purdue
    #
    #----   p.ncard: 'newcard'
    #       p.nvar: 6
    #       p.err: 1
    #       p.errmess: 'no errors detected on assignment statement'
    #       p.var: [1x6 struct]
    #          p.var(1).name: 'x'
    #          p.var(1).type: 'number'
    #          p.var(1).nval: 1
    #          p.var(1).val: 3
    #          p.var(2).name: 'rat'
    #          p.var(2).type: 'string'
    #          p.var(2).nval: 1
    #          p.var(2).val: {'cat'}
    #          p.var(3).name: 'alpha'
    #          p.var(3).type: 'number'
    #          p.var(3).nval: 1
    #          p.var(3).val: 1.0000e+012
    #          p.var(4).name: 'wl'
    #          p.var(4).type: 'number'
    #          p.var(4).nval: 5
    #          p.var(4)val: [0.3000 0.4000 0.5000 inf inf]
    #          p.var(5).name: 'string'
    #          p.var(5).type: 'string'
    #          p.var(5).nval: 4
    #          p.var(5).val: {'ab'  'cdcd'  'fgh'  'z'}
    #          p.var(6).name: 'purdue'
    #          p.var(6).type: 'empty'
    #          p.var(6).nval: 0
    #          p.var(6).val: []
    #
    #####   *mess  hello world!
    #
    #----   p.ncard: '*mess'
    #       p.err: 2
    #       p.errmess: 'no errors detected on title statement'
    #       p.title: 'hello world!  '
    #
    #####   test2 f=false r=3/yes
    #
    #----   p.ncard: 'test2'
    #       p.err: 999
    #       p.errmess: 'cannot mix numbers and strings in: r'
    #       p.nvar: 2
    #       p.var: [1x2 struct]
    #          p.var(1).name: 'f'
    #          p.var(1).type: 'string'
    #          p.var(1).nval: 1
    #          p.var(1).val: {'false'}
    #          p.var(2).name: 'r'
    #          p.var(2).type: 'number'
    #          p.var(2).nval: 2
    #          p.var(2).val: 3

    # Read first line of statement definition - ignore comment lines
    line = []
    line.append(' ')

    while line[0] == ' ':
        if fin.read=='':
            p.ncard = 'End-of-File'
            p.err = -1
            p.errmess = 'end of file'
            return

        line = next(fin, None)
        if line == None:
            p.err = -1
            p.errmess = 'return due to none'
            #p.ncard = 'End-of-File'
            #p.err = -1
            #p.errmess = 'end of file'
            return


        if line:
            if (np.size(list(line))== 1  or line[0] == '$' or line[0] == '%'):
                line_temp = list(line)
                line_temp[0] = ' '
                line = ''.join(line_temp)
        else:
            line = ' '

    statement = line

    # Read remaining lines of definition statement if continued

    line_temp = list(line)
    line_temp[0] = ' '
    line = ''.join(line_temp)
    while (line[0] == ' '):

        line = next(fin, None)
        if line == None:
            break

        if np.size(list(line))== 1:
            line_temp = list(line)
            line_temp[0] = ' '
            line = ''.join(line_temp)

        if line:
            if(line[0]=='+'):
                line_temp = list(line)
                line_temp[0] = ' '
                line = ''.join(line_temp)
                statement += line
            else:
                #fin.seek((np.size(line)+1), 1)
                break
        else:
           break
    #      line[0]=' '

    # replace all occurrences of ',' with a space, ' '

    statement = statement.replace(',', ' ')

    # decode the statement

    # check for title statement

    if statement[1] == '*':
        [p.ncard, rem] = statement.split(' ',1)
        p.title = rem.lstrip()
        p.err = 2
        p.errmess = 'no errors detected on title statement'
        return

    # get statement name
    if len(statement.split(' ',1)) != 1:
        [p.ncard, rem] = statement.split(' ', 1)
        rem = rem.rstrip()
    else:
        p.ncard = statement
        p.ncard = p.ncard.rstrip()
        rem = ''

    # If the end card is encounterd exit out from the parser
    if p.ncard == 'end':
        p.err = -1
        return

    # extract tokens

    i=0

    while i >= 0:
        if rem is None:
            break
        rem = rem.lstrip()

        try:
            [token, rem] = rem.split(' ',1)
        except ValueError:
            token = rem
            rem = None
        n = np.size(token)
        if(n==0):
           break

        i += 1
        if i == 1:
           list1 = [token]
        else:
           list1.append(token)

    # interpret list of tokens
    n = np.size(list1)
    p.nvar = n
    for ind in np.arange(n):
        p.add_var()

    for i in np.arange(0,n):
        [vname, rem] = list1[i].split('=',1)
        rem = rem.lstrip()
        rem = rem.rstrip()
        p.var[i].name = vname.rstrip()
        p.var[i].type='      '
        if rem:
            # test for array of values
            delims = rem.find('/')
            if delims == -1:
                nv = 1
            else:
                nd = np.size(delims)
                nv = nd+1
            p.var[i].nval = nv
            #p.var[i].add_val(nv)
            for j in np.arange(0,nv):
                if len(rem.split('/')) != 1:
                    [value, rem] = rem.split('/',1)
                else:
                    value = rem
                rem = rem.lstrip()
                rem = rem.rstrip()
                try:
                    number = float(value)
                except ValueError:
                    number = ''
                if number == '':
                    if p.var[i].type == 'number':
                        p.err = 999
                        p.errmess = 'cannot mix numbers and strings in:' + p.var[i].name
                        return
                    p.var[i].type = 'string'
                    p.var[i].val = value
                else:
                    if p.var[i].type == 'string':
                        p.err = 999
                        p.errmess = 'cannot mix numbers and strings in:' + p.var[i].name

                    p.var[i].type = 'number'
                    p.var[i].val = number
        else:
            p.var[i].nval = 0
            p.var[i].type = 'empty'

    p.err = 1
    p.errmess = 'no errors detected on assignment statement'

    return

