#!MC 1410
# Phasefield 3D (WBM): figura de SOLUTO (Cu).
#   - iso-superficie PhaseField = 0.5 (interface) colorida pela concentracao
#     de soluto (variavel 6, Conc_Solute), com iluminacao Gouraud; sem cortes
# Exporta PNG de 2000 px de largura.
#
# Uso (batch, Tecplot 360 EX 2017 R2 ou mais novo):
#   "C:\Program Files\Tecplot\Tecplot 360 EX 2017 R2\bin\tec360.exe" -b -p dendrita_wbm.mcr
# ou pela interface: Scripting > Play Macro/Script.
#
# >>> EDITE AS DUAS LINHAS ABAIXO <<<
#   ARQ = arquivo de saida do phasefield (output/saida_<IT>_time_<t>s.dat)
#   PNG = imagem a gerar
# O Tecplot do Windows nao le o sistema de arquivos do WSL: copie antes o arquivo
# para uma pasta do Windows (ex.: cp ~/run/output/saida_... /mnt/c/Users/<voce>/...).
$!VarSet |ARQ| = 'C:\pf_resultados\saida_wbm.dat'
$!VarSet |PNG| = 'C:\pf_resultados\saida_wbm_3d.png'

$!ReadDataSet  '"|ARQ|"'
  ReadDataOption = New
  ResetStyle = Yes
$!PlotType = Cartesian3D

$!FieldLayers ShowMesh = No
$!FieldLayers ShowContour = No
$!FieldLayers ShowShade = No
$!FieldLayers ShowEdge = No
$!FieldLayers ShowIsosurfaces = Yes

# grupo 1 define a iso-superficie (PhaseField), grupo 2 colore (Cu)
$!SetContourVar Var = 4 ContourGroup = 1 LevelInitMode = ResetToNice
$!SetContourVar Var = 6 ContourGroup = 2 LevelInitMode = ResetToNice
$!IsoSurfaceAttributes 1  DefinitionContourGroup = 1
$!IsoSurfaceAttributes 1  IsoSurfaceSelection = OneSpecificValue
$!IsoSurfaceAttributes 1  Isovalue1 = 0.5
$!IsoSurfaceAttributes 1  Contour { Show = Yes  FloodColoring = Group2 }
$!IsoSurfaceAttributes 1  Shade { UseLightingEffect = Yes }
$!IsoSurfaceAttributes 1  Effects { LightingEffect = Gouraud }

$!ThreeDView PSIAngle = 65 ThetaAngle = 215 AlphaAngle = 0
$!View Fit
$!FrameLayout ShowBorder = No

$!ExportSetup ExportFormat = PNG
$!ExportSetup ImageWidth = 2000
$!ExportSetup ExportFName = '|PNG|'
$!Export
$!Quit
