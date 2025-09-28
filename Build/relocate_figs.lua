-- relocate_figs.lua
-- Removes only the first leading "../" from a relative image path.
function Image(el)
  -- The source path is in the 'src' field of the image element
  local path = el.src
  
  -- The pattern matches only the first instance of '../'
  local cleaned_path = path:gsub("^(%.%.%/)","" ,1)
  
  -- Update the image element with the new path
  el.src = cleaned_path
  
  -- Return the modified element
  return el
end
