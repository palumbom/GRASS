# License

```@eval
using Markdown, GRASS
license_file = joinpath(pkgdir(GRASS), "LICENSE.md")
if !isfile(license_file)
    license_file = joinpath(pkgdir(GRASS), "...", "LICENSE.md")
end
Markdown.parse_file(license_file)
```
