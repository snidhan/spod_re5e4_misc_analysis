function [saveFile] = getSaveFileName(saveFolder,custom_comment,calcMethod,Re,m,Nr,Nz,r_farf,z1,z2,SIGMA)

if isnumeric(SIGMA)
    SIGMA_real   = num2str(real(SIGMA));
    SIGMA_imag   = num2str(imag(SIGMA));
else
    SIGMA_real   = SIGMA;
    SIGMA_imag   = SIGMA;
end

saveFile   = ['LESmeanflowJet2D_' calcMethod '_Re=' num2str(Re) '_m=' num2str(m) '_Nr=' num2str(Nr)  '_r=' num2str(r_farf) '_Nz=' num2str(Nz) '_' num2str(z1) '>z>' num2str(z2) '_ReSIGMA=' SIGMA_real '_ImSIGMA=' SIGMA_imag];

if length(custom_comment)>0
    saveFile   = [saveFolder '/' saveFile  '_' custom_comment '.mat'];
else
    saveFile   = [saveFolder '/' saveFile '.mat'];
end
