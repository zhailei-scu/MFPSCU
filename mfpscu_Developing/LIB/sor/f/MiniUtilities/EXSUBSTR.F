        subroutine extract_substr( String, mxCount,count,substr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11!
!  Purpose: to extract substring indicated ""
!  Input:   STRING a string
!  Ouput:   COUNT  the number of substring found in the string
!           SUBSTR the substring founded
!
! Auther:   Hou Qing, Inst. of Nucl. Sci.and Tech., Sichuan Union University
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       implicit none
       integer::i, index1, count, mxCount
       character *(*)::String
       character *(*)::substr
       dimension substr(mxCount)
       integer::curchar

       count=0
       index1=0
       i = 1
       substr = ""
       do while(.true.)
          if( iachar(String(i:i)) .eq. curchar .and. index1.gt.0 ) then
             count = count + 1
             substr(count)=string(index1:i-1)
             index1 = 0
          else
             if(String(i:i) .eq. '"' .or. String(i:i).eq."'"  )then
                index1 = i+1
                curchar = iachar(String(i:i))
             end if
          end if

          i=i+1
          if(i .gt. len_trim(string) .or. count .ge. mxCount) then
             exit
          end if
       enddo
       return
       end subroutine extract_substr
