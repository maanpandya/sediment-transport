using Plots

create_gifs = false
plot_separation_of_omega_equation_type = false
create_frequency_jump_up_plot = true
create_frequency_jump_down_plot = true

function cut_lists(lst1, lst2)
    # Ensure both lists are of the same length
    if length(lst1) != length(lst2)
        throw(ArgumentError("Both lists must have the same length"))
    end

    # Iterate through the list and check the condition
    for i in 1:length(lst1)-1
        if abs(lst1[i]) > 2 * abs(lst1[i+1])
            # Cut both lists from the next element to the end
            return lst1[i+1:end], lst2[i+1:end]
        end
    end
    
    # If the condition was never met, return the original lists
    return lst1, lst2
end

function find_drop_index(lst)
    for i in 1:(length(lst) - 1)  # Iterate until second to last element
        if lst[i] < lst[i + 1]   # Check if the next element is smaller
            return i              # Return the index of the drop
        end
    end
    return nothing  # Return nothing if no drop is found
end

function find_closest_index(lst, target)
    # Find the index of the element closest to the target value
    return argmin(abs.(lst .- target))
end

#α -> linear stiffness coefficient
#γ -> amplitude of periodic driving force
#β -> non-linear stiffness coefficient
#δ -> damping coefficient
#ω -> angular frequency of periodic driving force 
#a -> amplitude of the steady-state response 

α = 1.0
β = 0.04
δ = 0.1
γ = 1
a = collect(range(0, 10, 100))

a_plot_pos = [] #Contains the values of the steady-state response amplitude to be ploted for the positive ω formula
a_plot_neg = [] #Contains the values of the steady-state response amplitude to be ploted for the negative ω formula
ω_plot_pos = [] #Contains the values of the steady-state response driving angular frequency from the positive formula to be ploted
ω_plot_neg = [] #Contains the values of the steady-state response driving angular frequency from the negative formula to be ploted

#Calculate ω from analytical solution
for i in 1:length(a)

    discriminant = ((γ^2)/(a[i]^2)) - δ^2*α + ((δ^4)/(4)) - ((3*δ^2*β*a[i]^2)/(4))

    #Make sure ω^2 is not complex
    if discriminant > 0

        ω_1 = α - ((δ^2)/(2)) + ((3*β*a[i]^2)/(4)) + sqrt(discriminant) #Analytical solution 1
        ω_2 = α - ((δ^2)/(2)) + ((3*β*a[i]^2)/(4)) - sqrt(discriminant) #Analytical solution 2
    
        #Only take positive/physical values of ω^2
        if ω_1 >= 0
            push!(a_plot_pos, a[i])
            push!(ω_plot_pos, sqrt(ω_1))    
        end

        if ω_2 >= 0
             push!(a_plot_neg, a[i])
             push!(ω_plot_neg, sqrt(ω_2))    
        end
    end
end

#Remove infinity values from the obtained lists

indices_to_remove_neg = findall(x -> isinf(x), ω_plot_neg)
indices_to_remove_pos = findall(x -> isinf(x), ω_plot_pos)

ω_plot_neg = deleteat!(ω_plot_neg, indices_to_remove_neg)
a_plot_neg = deleteat!(a_plot_neg, indices_to_remove_neg)
ω_plot_pos = deleteat!(ω_plot_pos, indices_to_remove_pos)
a_plot_pos = deleteat!(a_plot_pos, indices_to_remove_pos)


#Sort the plot so that it is plotted with increasing and decreasing angular frequency depending on the formula used 
#sorted_indices = sortperm(ω_plot)

#Reorder both lists according to these indices
#ω_sorted = ω_plot[sorted_indices]
#a_sorted = a_plot[sorted_indices]

#println(ω_sorted)
#println(a_sorted)

if plot_separation_of_omega_equation_type
    plot(ω_plot_pos, a_plot_pos, xlim=(0, 3), ylim=(0, 8.5), label=false, color="red", title="Single Harmonic Duffing Frequency Response (positive formula)", size=(700, 600), titlepadding=100)
    scatter!(ω_plot_pos, a_plot_pos, color="orange", label="Positive")
    plot!(ω_plot_neg, a_plot_neg, color="dark green", label=false)
    scatter!(ω_plot_neg, a_plot_neg, color="green", label="Negative")
    # savefig("C:\\Users\\LENOVO\\Desktop\\duffing_single_response.png")
end

sorted_indices_neg = sortperm(ω_plot_neg) 
ω_sorted_neg = ω_plot_neg[sorted_indices_neg]
a_sorted_neg = a_plot_neg[sorted_indices_neg]

sorted_indices_pos = sortperm(ω_plot_pos, rev=true)
ω_sorted_pos = ω_plot_pos[sorted_indices_pos]
a_sorted_pos = a_plot_pos[sorted_indices_pos]

#Obtain frequency jump down behavior
differences = abs.(ω_sorted_pos .- ω_sorted_neg[end])
combined_condition = [(differences[i], a_sorted_pos[i]) for i in 1:length(ω_sorted_pos)]
sorted_combined = sort(combined_condition, by = x -> (x[1], x[2]))
jump_down_index = findfirst(x -> x == sorted_combined[1], combined_condition)

ω_down_pos = reverse(ω_sorted_pos[1:jump_down_index])
a_down_pos = reverse(a_sorted_pos[1:jump_down_index])

a_down_pos, ω_down_pos = cut_lists(a_down_pos, ω_down_pos)

ω_down = vcat(ω_sorted_neg, ω_down_pos)
a_down = vcat(a_sorted_neg, a_down_pos)

#Obtain the frequency jump up behavior
frequency_reversal_index = find_drop_index(ω_plot_pos)
ω_up_pos = ω_plot_pos[1:frequency_reversal_index]
a_up_pos = a_plot_pos[1:frequency_reversal_index]

frequency_reversal_index_2 = find_closest_index(ω_plot_neg, ω_up_pos[end])
ω_up_neg = reverse(ω_plot_neg[1:frequency_reversal_index_2])
a_up_neg = reverse(a_plot_neg[1:frequency_reversal_index_2])

ω_up = vcat(ω_up_pos, ω_up_neg)
a_up = vcat(a_up_pos, a_up_neg)

#Join the two behaviors
ω_plot = vcat(ω_down, ω_up)
a_plot = vcat(a_down, a_up)

if create_frequency_jump_up_plot
    plot(ω_up, a_up, xlim=(0, 3), ylim=(0, 8.5), label=false, color="red", title="Single Harmonic Duffing\n Frequency Response for Decreasing ω", size=(700, 600), titlepadding=100)
    scatter!(ω_up, a_up, color="dark red", label=false, xlabel="ω", ylabel="(u^2 + v^2)^1/2")
    # savefig("C:\\Users\\LENOVO\\Desktop\\duffing_single_response_jump_up.png")
end

if create_frequency_jump_down_plot
    plot(ω_down, a_down, xlim=(0, 3), ylim=(0, 8.5), label=false, color="green", title="Single Harmonic Duffing\n Frequency Response for Increasing ω", size=(700, 600), titlepadding=100)
    scatter!(ω_down, a_down, color="dark green", label=false, xlabel="ω", ylabel="(u^2 + v^2)^1/2")
    # savefig("C:\\Users\\LENOVO\\Desktop\\duffing_single_response_jump_down.png")
end


if create_gifs
    #Animate hysteresis behavior
    anim1 = @animate for i in 1:length(ω_down)
        scatter(ω_down[1:i], a_down[1:i],xlim=(0, 3), ylim=(0, 8.5), label=false, title="Single Harmonic Duffing\n Frequency Response for Increasing ω", xlabel="ω", ylabel="(u^2 + v^2)^1/2", size(800, 600))
    end
    anim2 = @animate for i in 1:length(ω_up)
        scatter(ω_up[1:i], a_up[1:i],xlim=(0, 3), ylim=(0, 8.5), label=false, title="Single Harmonic Duffing\n Frequency Response for Decreasing ω", xlabel="ω", ylabel="(u^2 + v^2)^1/2", size(800, 600))
    end

    #Save as a GIF
    # gif(anim1, "C:\\Users\\LENOVO\\Desktop\\frequency_response_duffing_single_harmonic1.gif", fps=10)
    # gif(anim2, "C:\\Users\\LENOVO\\Desktop\\frequency_response_duffing_single_harmonic2.gif", fps=10)
end