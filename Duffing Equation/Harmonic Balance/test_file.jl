for i in 1:2:5
    println("Iteration $i")
end

lst = collect(1:12)
print(lst)
lst2 = filter(x -> x % 2 != 0, collect(1:21))
print(lst2)