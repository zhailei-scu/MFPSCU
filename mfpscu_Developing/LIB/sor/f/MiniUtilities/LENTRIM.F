       function LENTRIM(str) result(length)
       implicit none
       !---Dummy Vars---
       character(*),intent(in)::str
       integer,intent(out)::length
       !---Local Vars---
       integer::I
       integer::tempLength
       !---Body---
       tempLength = len_trim(str)


       DO I = tempLength,1,-1
        if(ICHAR(str(I:I)) .GT. 31 .and. str(I:I) .ne. CHAR(32)
     -   .and. str(I:I) .ne. CHAR(34)) then
            exit
        end if
       END DO

      length = I

      return
      end function
