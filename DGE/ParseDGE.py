import tabula
import pandas as pd
import os
from os import listdir
import re
import requests
from datetime import datetime
from glob import glob
from bs4 import BeautifulSoup
from urllib.parse import urljoin

##################################################
#     DESCARGA DE ARCHIVOS DE LA DGE
##################################################

#Download files
#https://www.gob.mx/salud/documentos/informes-semanales-para-la-vigilancia-epidemiologica-de-viruela-simica-en-mexico
url = 'https://www.gob.mx/salud/documentos/informes-semanales-para-la-vigilancia-epidemiologica-de-viruela-simica-en-mexico'
reqs = requests.get(url)
soup = BeautifulSoup(reqs.text, 'html.parser')

download_folder = os.path.join("DGE","reportes_pdf")
save_folder     = os.path.join("DGE","casos_csv")

#Directorio
fnames = [os.path.join(download_folder, f) for f in listdir(download_folder)]
fnames = [name for name in fnames if name.find(".pdf") != -1]


for linkname in soup.find_all("a"):

    if (linkname.get("href") is not None) and (len(re.findall(r'[Ii]nforme_[Tt]ecnico_[Vv]iruela_[Ss]imica_.*.pdf', linkname.get('href'))) > 0):
        pdf = urljoin("https://www.gob.mx", linkname.get("href"))
        response = requests.get(pdf)
        fname = re.sub(r'https://.*/file/[0-9]+/[Ii]nforme_[Tt]ecnico_[Vv]iruela_[Ss]imica_|_|.pdf','',pdf)
        fname = os.path.join(download_folder, fname + ".pdf")
        # Write content in pdf file
        if (fname in fnames):
            print("File " + fname + " already downloaded")
        else:
            print("File " + fname + " writen")
            pdf = open(fname, 'wb')
            pdf.write(response.content)
            pdf.close()


##################################################
#     LECTURA DE ARCHIVOS DE PDF A CSV
##################################################

#Estados para verificar en el sistema
estados = {"CHIAPAS","VERACRUZ","NUEVO LEÓN","GUERRERO","MICHOACÁN","TAMAULIPAS",
           "JALISCO","MORELOS","PUEBLA","NAYARIT","TABASCO","QUINTANA ROO","SINALOA",
           "YUCATÁN","SONORA","MÉXICO","COLIMA","SAN LUIS POTOSÍ","OAXACA","HIDALGO",
           "COAHUILA","CAMPECHE","BAJA CALIFORNIA","BAJA CALIFORNIA SUR","CHIHUAHUA",
           "DURANGO","GUANAJUATO","QUERÉTARO","TLAXCALA","ZACATECAS","QUERETARO","SAN LUIS POTOSI",
           "AGUASCALIENTES","CIUDAD DE MÉXICO","BAJA CALIFORNIA **","ESTADO DE MÉXICO"}

#Directorio
fnames = [os.path.join(download_folder, f) for f in listdir(download_folder)]
fnames = [name for name in fnames if name.find(".pdf") != -1]



for name in fnames:
    print("Leyendo " + name)

    fecha =  datetime.strptime(re.sub(r'.*/|.pdf','',name), "%d%m%y")

    if fecha == datetime.strptime("150822", "%d%m%y"):
        coord = (210, 62, 680, 290)
    elif fecha == datetime.strptime("220822", "%d%m%y"):
        coord = (193, 65, 580, 290)
    elif fecha == datetime.strptime("290822", "%d%m%y"):
        coord = (225, 65, 610, 270)
    elif fecha == datetime.strptime("050922", "%d%m%y"):
        coord = (225, 65, 700, 270)
    elif fecha == datetime.strptime("041022", "%d%m%y"):
        coord = (159, 65, 700, 290)
    elif fecha == datetime.strptime("111022", "%d%m%y"):
        coord = (153, 68, 650, 378)
    elif fecha == datetime.strptime("270922", "%d%m%y"):
        coord = (184, 65, 680, 310)
    elif fecha == datetime.strptime("190922", "%d%m%y"):
        coord = (210, 62, 690, 304)
    elif fecha == datetime.strptime("080822", "%d%m%y"):
        coord = (237, 62, 470, 320)
    elif fecha == datetime.strptime("010822", "%d%m%y"):
        coord = (202, 62, 397, 320)
    elif fecha == datetime.strptime("120922", "%d%m%y"):
        coord = (195, 62, 680, 310)
    elif fecha == datetime.strptime("171022", "%d%m%y"):
        coord = (169, 62, 650, 378)
    else:
        coord = (169, 62, 650, 378)

    if fecha == datetime.strptime("251022", "%d%m%y"):
        d = {'Estado': ["Ciudad de México", "Yucatán", "Jalisco", "Quintana Roo", "Estado de México",
                        "Tabasco","Morelos","Querétaro","Campeche","Nayarit","Nuevo León","Baja California",
                        "Colima","Aguascalientes","Chiapas","Sinaloa","Hidalgo","Puebla","Tlaxcala","Chihuahua","Coahuila",
                        "Guanajuato","Tamaulipas","Veracruz","Oaxaca","Zacatecas","Sonora","Baja California Sur",
                        "Michoacán","Guerrero","San Luis Potosí","Durango"],
             'Casos': [1601, 93, 301, 58, 272, 39, 16, 17, 7, 9, 30, 19, 4, 7, 28, 13, 12, 25, 5, 13, 11, 16,
                        9, 20, 8, 3, 4, 1, 6, 4, 2, 1]
             }
        df = pd.DataFrame(d)
    else:
        table = tabula.read_pdf(name, pages=3, area=coord, pandas_options={'header': None})
        df = table[0]
        df = df[df[0].str.upper().isin(estados)]


    if (df.shape[1] == 2):
        df.columns = ["Estado","Casos"]
    elif (df.shape[1] == 3):
        df.columns = ["Estado","Casos","Tasa de incidencia"]
    elif (df.shape[1] == 6):
        df.columns = ["Estado", "Casos", "Tasa de incidencia","Ambulatorios","Hospitalizados","Altas"]
    df = df.assign(Fecha_reporte = fecha)

    if name != fnames[0]:
        df_all = pd.concat([df_all.reset_index(drop=True), df], axis = 0)
    else:
        df_all = df

df_all.to_csv(os.path.join(save_folder,"MPX_reporte_dge.csv"), index = False)
