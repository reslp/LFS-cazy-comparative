#!/bin/bash
set -e

echo "In 10 seconds this script will download assemblies from NCBI used for phylogenomic reconstructions. This can also be done with phylociraptor however the version of phylociraptor used in the paper does not support downloading specific accession. Hence, for reproducibility genome download is handled with this script."
sleep 10
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/814/355/GCA_022814355.1_ASM2281435v1/GCA_022814355.1_ASM2281435v1_genomic.fna.gz -O Acarospora_aff_strigata.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/814/335/GCA_022814335.1_ASM2281433v1/GCA_022814335.1_ASM2281433v1_genomic.fna.gz -O Agyrium_rufum.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/814/255/GCA_022814255.1_ASM2281425v1/GCA_022814255.1_ASM2281425v1_genomic.fna.gz -O Bachmanniomyces_sp_S44760.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/814/315/GCA_022814315.1_ASM2281431v1/GCA_022814315.1_ASM2281431v1_genomic.fna.gz -O Hypocenomyce_scalaris.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/814/295/GCA_022814295.1_ASM2281429v1/GCA_022814295.1_ASM2281429v1_genomic.fna.gz -O Icmadophila_ericetorum.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/814/265/GCA_022814265.1_ASM2281426v1/GCA_022814265.1_ASM2281426v1_genomic.fna.gz -O Lambiella_insularis.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/814/235/GCA_022814235.1_ASM2281423v1/GCA_022814235.1_ASM2281423v1_genomic.fna.gz -O Lignoscripta_atroalba.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/814/215/GCA_022814215.1_ASM2281421v1/GCA_022814215.1_ASM2281421v1_genomic.fna.gz -O Lobaria_immixta.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/814/175/GCA_022814175.1_ASM2281417v1/GCA_022814175.1_ASM2281417v1_genomic.fna.gz -O Loxospora_ochrophaea.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/814/195/GCA_022814195.1_ASM2281419v1/GCA_022814195.1_ASM2281419v1_genomic.fna.gz -O Mycoblastus_sanguinarius.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/814/155/GCA_022814155.1_ASM2281415v1/GCA_022814155.1_ASM2281415v1_genomic.fna.gz -O Peltigera_leucophlebia.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/814/125/GCA_022814125.1_ASM2281412v1/GCA_022814125.1_ASM2281412v1_genomic.fna.gz -O Pseudocyphellaria_aurata.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/814/105/GCA_022814105.1_ASM2281410v1/GCA_022814105.1_ASM2281410v1_genomic.fna.gz -O Ptychographa_xylographoides.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/814/085/GCA_022814085.1_ASM2281408v1/GCA_022814085.1_ASM2281408v1_genomic.fna.gz -O Puttea_exsequens.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/814/065/GCA_022814065.1_ASM2281406v1/GCA_022814065.1_ASM2281406v1_genomic.fna.gz -O Schaereria_dolodes.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/814/045/GCA_022814045.1_ASM2281404v1/GCA_022814045.1_ASM2281404v1_genomic.fna.gz -O Sticta_canariensis.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/814/025/GCA_022814025.1_ASM2281402v1/GCA_022814025.1_ASM2281402v1_genomic.fna.gz -O Stictis_urceolatum.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/813/995/GCA_022813995.1_ASM2281399v1/GCA_022813995.1_ASM2281399v1_genomic.fna.gz -O Thelotrema_lepadinum.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/813/965/GCA_022813965.1_ASM2281396v1/GCA_022813965.1_ASM2281396v1_genomic.fna.gz -O Toensbergia_leucococca.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/813/945/GCA_022813945.1_ASM2281394v1/GCA_022813945.1_ASM2281394v1_genomic.fna.gz -O Trapelia_coarctata.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/813/925/GCA_022813925.1_ASM2281392v1/GCA_022813925.1_ASM2281392v1_genomic.fna.gz -O Varicellaria_rhodocarpa.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/813/905/GCA_022813905.1_ASM2281390v1/GCA_022813905.1_ASM2281390v1_genomic.fna.gz -O Xylographa_bjoerkii.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/813/885/GCA_022813885.1_ASM2281388v1/GCA_022813885.1_ASM2281388v1_genomic.fna.gz -O Xylographa_carneopallida.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/813/865/GCA_022813865.1_ASM2281386v1/GCA_022813865.1_ASM2281386v1_genomic.fna.gz -O Xylographa_opegraphella.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/813/845/GCA_022813845.1_ASM2281384v1/GCA_022813845.1_ASM2281384v1_genomic.fna.gz -O Xylographa_pallens.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/813/825/GCA_022813825.1_ASM2281382v1/GCA_022813825.1_ASM2281382v1_genomic.fna.gz -O Xylographa_parallela.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/813/805/GCA_022813805.1_ASM2281380v1/GCA_022813805.1_ASM2281380v1_genomic.fna.gz -O Xylographa_soralifera.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/813/785/GCA_022813785.1_ASM2281378v1/GCA_022813785.1_ASM2281378v1_genomic.fna.gz -O Xylographa_trunciseda.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/813/765/GCA_022813765.1_ASM2281376v1/GCA_022813765.1_ASM2281376v1_genomic.fna.gz -O Xylographa_vitiligo.fna.gz

## previously published:

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/989/075/GCA_002989075.1_artRad1.0/GCA_002989075.1_artRad1.0_genomic.fna.gz -O Arthonia_radiata.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/336/255/GCA_003336255.1_ASM333625v1/GCA_003336255.1_ASM333625v1_genomic.fna.gz -O Aureobasidium_pullulans.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/143/535/GCA_000143535.4_ASM14353v4/GCA_000143535.4_ASM14353v4_genomic.fna.gz -O Botrytis_cinerea.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/585/585/GCA_000585585.1_Capr_coro_CBS_617_96_V1/GCA_000585585.1_Capr_coro_CBS_617_96_V1_genomic.fna.gz -O Capronia_coronata.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/978/885/GCA_000978885.1_ASM97888v1/GCA_000978885.1_ASM97888v1_genomic.fna.gz -O Ceratocystis_platani.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/521/265/GCA_003521265.1_NYBG_CetLin_1.0/GCA_003521265.1_NYBG_CetLin_1.0_genomic.fna.gz -O Cetradonia_linearis.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/444/155/GCA_000444155.1_Clmac_v1/GCA_000444155.1_Clmac_v1_genomic.fna.gz -O Cladonia_macilenta.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/482/085/GCA_000482085.2_C_metacorallifera_KoLRI002260_v2/GCA_000482085.2_C_metacorallifera_KoLRI002260_v2_genomic.fna.gz -O Cladonia_metacorallifera.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/927/785/GCA_002927785.1_ASM292778v1/GCA_002927785.1_ASM292778v1_genomic.fna.gz -O Cladonia_uncialis.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/365/165/GCA_000365165.2_Clad_carr_CBS_160_54_V1/GCA_000365165.2_Clad_carr_CBS_160_54_V1_genomic.fna.gz -O Cladophialophora_carrionii.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/149/035/GCA_000149035.1_C_graminicola_M1_001_V1/GCA_000149035.1_C_graminicola_M1_001_V1_genomic.fna.gz -O Colletotrichum_graminicola.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/019/895/GCA_003019895.1_Coniella_lustricola_v1.0/GCA_003019895.1_Coniella_lustricola_v1.0_genomic.fna.gz -O Coniella_lustricola.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/879/275/GCA_001879275.1_Conlig_v1.0/GCA_001879275.1_Conlig_v1.0_genomic.fna.gz -O Coniochaeta_ligniaria.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/618/795/GCA_900618795.1_Astra/GCA_900618795.1_Astra_genomic.fna.gz -O Cyanodermella_asteris.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/800/745/GCA_000800745.1_ASM80074v1/GCA_000800745.1_ASM80074v1_genomic.fna.gz -O Diaporthe_longicolla.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/240/705/GCA_002240705.1_PX439/GCA_002240705.1_PX439_genomic.fna.gz -O Elaphomyces_granulatus.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/464/535/GCA_000464535.1_EPUS/GCA_000464535.1_EPUS_genomic.fna.gz -O Endocarpon_pusillum.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/308/955/GCA_000308955.1_ETyphina_1.0/GCA_000308955.1_ETyphina_1.0_genomic.fna.gz -O Epichloe_typhina.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/798/715/GCA_000798715.1_ASM79871v1/GCA_000798715.1_ASM79871v1_genomic.fna.gz -O Erysiphe_necator.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/184/365/GCA_003184365.1_ASM318436v1/GCA_003184365.1_ASM318436v1_genomic.fna.gz -O Evernia_prunastri.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/651/985/GCA_001651985.1_ASM165198v1/GCA_001651985.1_ASM165198v1_genomic.fna.gz -O Fonsecaea_erecta.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/149/955/GCA_000149955.2_ASM14995v2/GCA_000149955.2_ASM14995v2_genomic.fna.gz -O Fusarium_oxysporum.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/337/335/GCA_905337335.1_GTX0158_lecanoromycete/GCA_905337335.1_GTX0158_lecanoromycete_genomic.fna.gz -O Gomphillus_americanus.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/442/125/GCA_000442125.1_Cafla_v1/GCA_000442125.1_Cafla_v1_genomic.fna.gz -O Gyalolechia_flavorubescens.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/573/585/GCA_002573585.1_Heli_gris_5409_V1/GCA_002573585.1_Heli_gris_5409_V1_genomic.fna.gz -O Helicocarpus_griseus.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/775/035/GCA_002775035.1_ASM277503v1/GCA_002775035.1_ASM277503v1_genomic.fna.gz -O Hypoxylon_pulicicidum.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/319/055/GCA_002319055.1_ASM231905v1/GCA_002319055.1_ASM231905v1_genomic.fna.gz -O Knufia_petricola.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/254/425/GCA_003254425.1_ASM325442v1/GCA_003254425.1_ASM325442v1_genomic.fna.gz -O Lasallia_hispanica.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/169/345/GCA_900169345.1_Lasallia_pustulata_v1/GCA_900169345.1_Lasallia_pustulata_v1_genomic.fna.gz -O Lasallia_pustulata.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/871/045/GCA_002871045.1_ASM287104v1/GCA_002871045.1_ASM287104v1_genomic.fna.gz -O Magnaporthe_grisea.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/128/795/GCA_900128795.2_FCH_10_5/GCA_900128795.2_FCH_10_5_genomic.fna.gz -O Malbranchea_cinnamomea.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/184/285/GCA_003184285.1_ASM318428v1/GCA_003184285.1_ASM318428v1_genomic.fna.gz -O Monascus_purpureus.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/812/245/GCA_000812245.1_ASM81224v1/GCA_000812245.1_ASM81224v1_genomic.fna.gz -O Onygena_corvina.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/710/275/GCA_000710275.1_ASM71027v1/GCA_000710275.1_ASM71027v1_genomic.fna.gz -O Penicillium_chrysogenum.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/006/345/GCA_001006345.1_UCRPA7V1.0/GCA_001006345.1_UCRPA7V1.0_genomic.fna.gz -O Phaeomoniella_chlamydospora.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/299/255/GCA_001299255.1_ASM129925v1/GCA_001299255.1_ASM129925v1_genomic.fna.gz -O Phialophora_attae.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/604/955/GCA_001604955.1_ASM160495v1/GCA_001604955.1_ASM160495v1_genomic.fna.gz -O Phyllosticta_citricarpa.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/226/545/GCA_000226545.1_ASM22654v1/GCA_000226545.1_ASM22654v1_genomic.fna.gz -O Podospora_anserina.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/184/345/GCA_003184345.1_ASM318434v1/GCA_003184345.1_ASM318434v1_genomic.fna.gz -O Pseudevernia_furfuracea.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/868/215/GCA_003868215.1_ASM386821v1/GCA_003868215.1_ASM386821v1_genomic.fna.gz -O Pseudophaeomoniella_oleicola.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/073/195/GCA_003073195.1_RamPxa02_v1.0/GCA_003073195.1_RamPxa02_v1.0_genomic.fna.gz -O Ramalina_intermedia.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/956/345/GCA_001956345.1_RamPxa01_v1/GCA_001956345.1_RamPxa01_v1_genomic.fna.gz -O Ramalina_peruviana.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/182/805/GCA_000182805.2_ASM18280v2/GCA_000182805.2_ASM18280v2_genomic.fna.gz -O Sordaria_macrospora.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/829/775/GCA_000829775.1_TcelY94_1.0/GCA_000829775.1_TcelY94_1.0_genomic.fna.gz -O Talaromyces_cellulolyticus.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/599/835/GCA_001599835.1_JCM_12817_assembly_v001/GCA_001599835.1_JCM_12817_assembly_v001_genomic.fna.gz -O Thermoascus_crustaceus.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/513/885/GCA_001513885.1_ASM151388v1/GCA_001513885.1_ASM151388v1_genomic.fna.gz -O Thielaviopsis_musarum.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/167/675/GCA_000167675.2_v2.0/GCA_000167675.2_v2.0_genomic.fna.gz -O Trichoderma_reesei.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/611/775/GCA_000611775.1_Umbmu1.0/GCA_000611775.1_Umbmu1.0_genomic.fna.gz -O Umbilicaria_muehlenbergii.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/003/515/GCA_000003515.2_ASM351v2/GCA_000003515.2_ASM351v2_genomic.fna.gz -O Uncinocarpus_reesii.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/150/675/GCA_000150675.2_ASM15067v2/GCA_000150675.2_ASM15067v2_genomic.fna.gz -O Verticillium_dahliae.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/353/285/GCA_004353285.1_ASM435328v1/GCA_004353285.1_ASM435328v1_genomic.fna.gz -O Xylaria_grammica.fna.gz

echo
echo "Moving and unpacking genomes"
mv *.gz assemblies
gunzip assemblies/*.gz

echo
echo
echo "IMPORTANT:"
echo "The Cladonia grayi genome is not available on NCBI. Please download the Cladonia grayi genome manually, you have to create a free account on the JGI website to do so:"
echo "https://mycocosm.jgi.doe.gov/Clagr3/Clagr3.info.html"
echo "All other genomes have been downloaded and unpacked."







