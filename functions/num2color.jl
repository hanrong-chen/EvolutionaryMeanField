################################################################################
# 6/2/16: extended color scheme
# 9/9/16: make colordict an array
# 9/9/16: remove colors from consideration (so that multiplicity is different)
# 9/9/16: handle num<=0
# 9/9/16: added colors
################################################################################
function num2color(num; range=30)
    colordict0=["blue", "green", "red", "cyan", "forestgreen",
    			"yellow", "darkblue", "lightgreen", "maroon", "purple",
                "tan", "gold", "plum", "olive", "skyblue",
                "teal", "brown", "turquoise", "crimson", "darkgreen",
                "khaki", "orange", "grey", "violet", "magenta",
                "limegreen", "silver", "salmon", "royalblue", "black"]
    
    # 9/9/16: colors that are 
    # too light (5): pink, lightblue, lavender, aqua, beige
    # not available (6): light purple, mauve, dark purple, bright green, navy blue, lilac

    # 9/9/16: remove colors from consideration (so that multiplicity is different)
    colordict=colordict0[1:range]
    
    # 9/9/16: handle num<=0
    while num<=0
        num+=range
    end

    return colordict[(num-1)%range+1]
end

# for i=1:20
#     vlines(i,0,1,color=num2color(i))
# end


################################################################################
# 9/1/16: num2color using black-box function from Colors package
# note 1: I'm getting rid of white as one of the returned colors
# note 2: num%num_colors==0 should return blue
################################################################################
# using Colors
# using Images

# num_colors=25
# initial_colors=[colorant"white",colorant"blue",colorant"green",colorant"red",
#               colorant"cyan",colorant"magenta",colorant"yellow",colorant"black"]
# # get rid of white
# color_list=separate(distinguishable_colors(num_colors+1,initial_colors))[2:end,:]

# function num2color(num)
#   this_num=num%num_colors+1
#   return [color_list[this_num,1][1],color_list[this_num,2][1],
#           color_list[this_num,3][1]]
# end