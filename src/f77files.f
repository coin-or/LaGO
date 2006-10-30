*     Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
*     All Rights Reserved.
*     This code is published under the Common Public License.
*
*     Author: Stefan Vigerske

      subroutine f77openfile(id, string, len)

      integer       id
      character*1   string(len)
      integer       len

      character*256 scopy

      scopy(1:256)=' '

      do i=1, len
        scopy(i:i)=string(i)
      enddo

      OPEN(id, FILE=scopy)

      end

*     ----------------------------------------------------------

      subroutine f77closefile(id)

      CLOSE(id)

      end
