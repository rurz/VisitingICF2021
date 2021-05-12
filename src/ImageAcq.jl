"This script is intended to manipulate the multidimensional arrays from the image transformation"

import Images

export imgload, chsplit

"`imgload(path_string)` returns the array of values of a monochromatic or polychromatic image. In the first case is onlye an N×N matrix, in the second case is a collection of three N×N matrices."
function imgload(path::String)
    return Images.load(path)
end

"`chsplit(img, c)` returns the c-channel of and image _img_. It requires a polychromatic imported image from the `imgload` function, and a channel c ∈ {1, 2, 3} whose equivalences are 1 = r (red), 2 = g (green), 3 = b (blue)."
function chsplit(img, c)
    chs = Images.channelview(img)
    if c == 1
        return chs[1,:,:]
    elseif c == 2
        return chs[2,:,:]
    elseif c == 3
        return chs[3,:,:]
    end
end

"`imgrep(r, g, b)` takes the channel arrays (original or transformated), and returns a stacked representation of dimensions N×N×3, suitable for plotting in the `imshow` command of PyPlot"
function imgrep(r,g,b)
    return permutedims(Images.StackedView(r, g, b), [3, 2, 1])
end
