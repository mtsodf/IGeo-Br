-- Filtro pandoc para os relatórios em PDF (usado por docs/relatorio.mk).
--  * troca as figuras .svg (usadas no GitHub) pelas versões .pdf geradas
--    por gerar_figuras.py, que o LaTeX inclui sem conversão;
--  * faz as imagens ocuparem a largura do texto;
--  * junta cada imagem com a sua legenda numa figura do LaTeX;
--  * links relativos para arquivos do repositório (../../tests/...) não
--    funcionam fora do GitHub: mantém só o texto do link.
function Image(img)
  img.src = img.src:gsub("%.svg$", ".pdf")
  img.attributes.width = img.attributes.width or "100%"
  return img
end

function Link(link)
  if not link.target:match("^%a[%w+.-]*:") and not link.target:match("^#") then
    return link.content
  end
end

-- A seta → não existe na maioria das fontes de texto: usa a do modo matemático.
function Str(s)
  if not s.text:find("→", 1, true) then return nil end
  local out, resto = {}, s.text
  while true do
    local i, j = resto:find("→", 1, true)
    if not i then break end
    if i > 1 then out[#out + 1] = pandoc.Str(resto:sub(1, i - 1)) end
    out[#out + 1] = pandoc.Math("InlineMath", "\\rightarrow")
    resto = resto:sub(j + 1)
  end
  if #resto > 0 then out[#out + 1] = pandoc.Str(resto) end
  return out
end

-- No README, cada imagem é seguida da legenda num parágrafo em itálico
-- ("*Figura N — ...*"). Junta os dois numa figura do LaTeX, para que imagem e
-- legenda não fiquem em páginas diferentes. O número já está no texto da
-- legenda; o rótulo automático "Figura N:" é desligado em metadados.yaml.
local function so_imagem(b)
  return b and b.t == "Para" and #b.content == 1 and b.content[1].t == "Image"
end

local function legenda(b)
  return b and b.t == "Para" and #b.content == 1 and b.content[1].t == "Emph"
end

function Blocks(blocos)
  local out = pandoc.Blocks{}
  local i = 1
  while i <= #blocos do
    if so_imagem(blocos[i]) and legenda(blocos[i + 1]) then
      out:insert(pandoc.Figure(pandoc.Plain{blocos[i].content[1]},
                               {long = {pandoc.Plain(blocos[i + 1].content[1].content)}}))
      i = i + 2
    else
      out:insert(blocos[i])
      i = i + 1
    end
  end
  return out
end

-- Os separadores "---" entre as seções do README não são necessários no PDF.
function HorizontalRule()
  return {}
end
