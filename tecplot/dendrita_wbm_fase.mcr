#!MC 1410
# Phasefield 3D (WBM): figura de FASE.
#   - iso-superficie PhaseField = 0.5 (interface) em vermelho solido, com iluminacao Gouraud
#   - cortes (slices) nas paredes x = 0 e y = 0 coloridos pelo PhaseField
#     (azul = liquido, vermelho = solido)
# Exporta PNG de 2000 px de largura.
#
# Uso (batch, Tecplot 360 EX 2017 R2 ou mais novo):
#   "C:\Program Files\Tecplot\Tecplot 360 EX 2017 R2\bin\tec360.exe" -b -p dendrita_wbm_fase.mcr
# ou pela interface: Scripting > Play Macro/Script.
#
# >>> EDITE AS DUAS LINHAS ABAIXO <<<
#   ARQ = arquivo de saida do phasefield (output/saida_<IT>_time_<t>s.dat)
#   PNG = imagem a gerar
# O Tecplot do Windows nao le o sistema de arquivos do WSL: copie antes o arquivo
# para uma pasta do Windows (ex.: cp ~/run/output/saida_... /mnt/c/Users/<voce>/...).
$!VarSet |ARQ| = 'C:\pf_resultados\saida_wbm.dat'
$!VarSet |PNG| = 'C:\pf_resultados\saida_wbm_fase_3d.png'

$!ReadDataSet  '"|ARQ|"'
  ReadDataOption = New
  ResetStyle = Yes
$!PlotType = Cartesian3D

$!FieldLayers ShowMesh = No
$!FieldLayers ShowContour = No
$!FieldLayers ShowShade = No
$!FieldLayers ShowEdge = No
$!FieldLayers ShowIsosurfaces = Yes

# grupo 1 = PhaseField, niveis de 0 a 1
$!SetContourVar Var = 4 ContourGroup = 1 LevelInitMode = ResetToNice
$!ContourLevels New ContourGroup = 1 RawData
21
0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1

# iso-superficie da interface, vermelha (solido), com iluminacao
$!IsoSurfaceAttributes 1  DefinitionContourGroup = 1
$!IsoSurfaceAttributes 1  IsoSurfaceSelection = OneSpecificValue
$!IsoSurfaceAttributes 1  Isovalue1 = 0.5
$!IsoSurfaceAttributes 1  Contour { Show = No }
$!IsoSurfaceAttributes 1  Shade { Show = Yes  Color = Red  UseLightingEffect = Yes }
$!IsoSurfaceAttributes 1  Effects { LightingEffect = Gouraud }

# cortes nas paredes do fundo coloridos pelo PhaseField
$!SliceLayers Show = Yes
$!SliceAttributes 1  SliceSurface = XPlanes
$!SliceAttributes 1  PrimaryPosition { X = 1E-09 }
$!SliceAttributes 1  Contour { Show = Yes  ContourType = Flood  FloodColoring = Group1 }
$!SliceAttributes 2  Show = Yes
$!SliceAttributes 2  SliceSurface = YPlanes
$!SliceAttributes 2  PrimaryPosition { Y = 1E-09 }
$!SliceAttributes 2  Contour { Show = Yes  ContourType = Flood  FloodColoring = Group1 }

$!ThreeDView PSIAngle = 65 ThetaAngle = 215 AlphaAngle = 0
$!View Fit
$!FrameLayout ShowBorder = No

$!ExportSetup ExportFormat = PNG
$!ExportSetup ImageWidth = 2000
$!ExportSetup ExportFName = '|PNG|'
$!Export
$!Quit
