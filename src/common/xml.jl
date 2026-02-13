"""
Some UPF files seem to have unescaped ampersands which cause XML parsing to fail.

This could also be a problem for other XML-based formats (e.g. PSML).

In XML, `&` actually _begins_ escape sequences (e.g. `<` is `&lt;`),
so we only want to replace `&`s that are not followed by a valid escape sequence.
"""
escape_ampersand(text) = replace(text, r"&(?!(amp|apos|quot|lt|gt);)" => "&amp;")
