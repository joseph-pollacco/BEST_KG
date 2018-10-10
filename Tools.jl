module array
   export SEARCH_INDEX, SEARCH_INDEX_STRING


   function SEARCH_INDEX_STRING(Array, SearchValue)
      N = length(Array)
      iSearchValue = 1
      
      for i in 1:N            
         if Array[i] == SearchValue
            return iSearchValue = i
            break
         end
      end
   end 



   function SEARCH_INDEX(Array, SearchValue)
      N = length(Array)
      iSearchValue = 1
      Value_SearchValue=1.
      Err_2 = 100000000000000.
      
      for i in 1:N
         Err_1 = abs(Array[i] - SearchValue)
         
         if Err_1 < Err_2
            iSearchValue = i
            Value_SearchValue = Array[i]
            Err_2 = Err_1
         end
      end
      return iSearchValue
   end 
end

