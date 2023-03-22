function TEB = TEB(t0,signal,Ns,info_binaire)
    signal_ech=signal(t0:Ns:end);
    signes_sorties=sign(signal_ech);
    info_bin_sortie=(signes_sorties+1)/2;
    vect_erreur=(info_binaire~=info_bin_sortie);
    TEB=mean(vect_erreur);
    
end