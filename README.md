# tmbFuns

A repository for my TMB function libraries. Part of my push to refactor my workflow to a more efficient state.

I hope these functions are useful beyond my own work. Feel free to suggest functions, or to fork and golf my code. Submit a pull request, I'd love to see this grow from community involvement.

The main thing I ask if you're golfing my code is that you refrain from changing function inputs and outputs (including side effects - I know, poor functional style), as this will create errors in my work where I apply these functions. Of course, I'm open to suggestion and discussion on all parts of this repo, so if you feel strongly about changing the function I/O start an issue.

## How to use

Add an #include statement to the top of the *.cpp file containing your TMB model. For example, to include the stockAssessmentFuns.hpp library:
    
    #include "<path-to>/tmbFuns/stockAssessmentFuns.hpp"

where <path-to> is replaced with the actual file path. 

Alternatively, you can copy your library file to your working directory, and use the statement
    
    #include "stockAssessmentFuns.hpp"

I don't recommend the second approach, though, as it will not update when or if you update the repository version of the library. This will potentially lead to multiple conflicting version across your system, which is the exact kind of behaviour this repository was created to avoid. As long as the function inputs and outputs for a function remain the same, a centralised library is a better option.

### Adding tmbFuns/ to compiler's include path

There is a way to add tmbFuns to the default include path for your compiler. If you use GCC, there is a stack-exchange post [here.](https://stackoverflow.com/questions/558803/how-to-add-a-default-include-path-for-gcc-in-linux) Although it looks like it's as simple as updating the paths in your base profile, I haven't attempted this. The reason is that I control my TMB models from an R script, so I anticipate some fussing with how R loads the bash environment.
