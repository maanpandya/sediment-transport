function find_drop_index(lst)
    for i in 1:(length(lst) - 1)  # Iterate until second to last element
        if lst[i] > lst[i + 1]   # Check if the next element is smaller
            return i              # Return the index of the drop
        end
    end
    return nothing  # Return nothing if no drop is found
end

# Example usage
lst = [1, 3, 5, 4, 6, 8]
index = find_drop_index(lst)
println(index)  # Output: 3 (where 5 > 4)