import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import TLorentzVector
import math
from math import *
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.wmass.recoZproducer import format_name

def evaluateMt(mu_pt, mu_phi, mu_mass, met_pt, met_phi) :
    muvec = TLorentzVector()
    muvec.SetPtEtaPhiM(mu_pt,0.0,mu_phi, mu_mass)
    metvec = TLorentzVector()
    metvec.SetPtEtaPhiM(met_pt,0.0,met_phi,0.0)
    hvec = -(muvec + metvec)
    mt, h_pt, h_phi  = hvec.M(), hvec.Pt(), hvec.Phi()
    met_par, met_per = metvec.Vect().Dot(muvec.Vect().Unit()), metvec.Vect().Cross(muvec.Vect().Unit()).Pz()
    h_par, h_per = hvec.Vect().Dot(muvec.Vect().Unit()), hvec.Vect().Cross(muvec.Vect().Unit()).Pz()
    return  (mt, h_pt, h_phi, met_par, met_per,h_par,h_per)

class additionalVariables(Module):
    def __init__(self, isMC=True, mudict={}, metdict={}, lightOutput=0):
        self.mudict = mudict
        self.metdict = metdict
        self.lightOutput = lightOutput
        self.isMC = isMC
        pass

    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for key_mu,mu in self.mudict.items() :
            for ivar_mu,var_mu in enumerate(mu["systs"]):
                for key_met,met in self.metdict.items() :
                    for ivar_met,var_met in enumerate(met["systs"]):
                        # skip double-systematic variation
                        if ("Up" in var_mu or "Down" in var_mu) and ("Up" in var_met or "Down" in var_met):
                            continue
                        # skip mixed reco-gen objects
                        if (key_met=="GEN" and key_mu!="GEN") or (key_met!="GEN" and key_mu=="GEN"):
                            continue
                        title_label = "MET: "+met["tag"]+" (syst: "+var_met+"), Muon: "+mu["tag"]+" (syst: "+var_mu+")"
                        branch_label_mu = format_name(mu["tag"], var_mu, met["tag"], var_met, False)
                        branch_label_Z  = ("GenFromRecoZ" if key_mu=="GEN" else "RecoZ")
                        branch_label_Z += format_name(mu["tag"], var_mu, met["tag"], var_met, True)
                        self.out.branch("%s_mt"   % branch_label_mu, "F", lenVar="nMuons", title="Mt, "+title_label)
                        if self.lightOutput<2:
                            self.out.branch("%s_hpt"  % branch_label_mu, "F", lenVar="nMuons", title="H pt, "+title_label)
                            self.out.branch("%s_hphi" % branch_label_mu, "F", lenVar="nMuons", title="H phi, "+title_label)
                            self.out.branch("%s_Wlikemt" % branch_label_mu, "F",  lenVar="nMuons", title="W-like Mt, "+title_label)
                            self.out.branch("%s_Wlikehpt" % branch_label_mu, "F",  lenVar="nMuons", title="W-like H pt, "+title_label)
                            self.out.branch("%s_Wlikehphi" % branch_label_mu, "F",  lenVar="nMuons", title="W-like H phi, "+title_label)
                            self.out.branch("%s_WlikeuPar" % branch_label_mu, "F",  lenVar="nMuons", title="W-like MET || mu, "+title_label)
                            self.out.branch("%s_WlikeuPer" % branch_label_mu, "F",  lenVar="nMuons", title="W-like MET T mu, "+title_label)
                            self.out.branch("%s_metPt" % branch_label_mu, "F", lenVar="nMuons", title="MET Pt, "+title_label)
                            self.out.branch("%s_uPar" % branch_label_mu, "F", lenVar="nMuons", title="MET || mu, "+title_label)
                            self.out.branch("%s_uPer" % branch_label_mu, "F", lenVar="nMuons", title="MET T mu, "+title_label)
                            self.out.branch("%s_hPar" % branch_label_mu, "F", lenVar="nMuons", title="H pt || mu, "+title_label)
                            self.out.branch("%s_hPer" % branch_label_mu, "F", lenVar="nMuons", title="H pt T mu, "+title_label)
                            self.out.branch("%s_WlikehPar" % branch_label_mu, "F", lenVar="nMuons", title="W-like H pt || mu, "+title_label)
                            self.out.branch("%s_WlikehPer" % branch_label_mu, "F", lenVar="nMuons", title="W-like H pt T mu, "+title_label)
                            self.out.branch("%s_uPar" % branch_label_Z, "F", title=branch_label_Z+" MET || Z, "+title_label)
                            self.out.branch("%s_uPer" % branch_label_Z, "F", title=branch_label_Z+" MET T Z, "+title_label)
                            self.out.branch("%s_hPar" % branch_label_Z, "F", title=branch_label_Z+" H pt || Z, "+title_label)
                            self.out.branch("%s_hPer" % branch_label_Z, "F", title=branch_label_Z+" H pt T Z, "+title_label)

        #self.out.branch("Muon_pfIso04", "F", lenVar="nMuons", title="muon absolute isolation")

        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        #collections
        muons = Collection(event, "Muon")
        if self.isMC:
            genparticles = Collection(event, "GenPart")

        #final vectors
        Mt_vec         = [0.0]*len(muons)
        Recoil_pt_vec  = [0.0]*len(muons)
        Recoil_phi_vec = [0.0]*len(muons)
        met_par_vec    = [0.0]*len(muons)
        met_per_vec    = [0.0]*len(muons)
        Recoil_par_vec    = [0.0]*len(muons)
        Recoil_per_vec    = [0.0]*len(muons)
        met_pt_vec = [0.0]*len(muons)
        Mt_Wlike_vec         = [0.0]*len(muons)
        Recoil_pt_Wlike_vec  = [0.0]*len(muons)
        Recoil_phi_Wlike_vec = [0.0]*len(muons)
        Recoil_par_Wlike_vec    = [0.0]*len(muons)
        Recoil_per_Wlike_vec    = [0.0]*len(muons)
        met_par_Wlike_vec    = [0.0]*len(muons)
        met_per_Wlike_vec    = [0.0]*len(muons)
        pfIso04_vec    = [0.0]*len(muons)

        for key_mu,mu in self.mudict.items() :
            for ivar_mu,var_mu in enumerate(mu["systs"]):
                for key_met,met in self.metdict.items() :
                    for ivar_met,var_met in enumerate(met["systs"]):
                        if ("Up" in var_mu or "Down" in var_mu) and ("Up" in var_met or "Down" in var_met):
                            continue
                        if (key_met=="GEN" and key_mu!="GEN") or (key_met!="GEN" and key_mu=="GEN"):
                            continue
                        branch_label_mu = format_name(mu["tag"], var_mu, met["tag"], var_met, False)
                        met_pt  = getattr(event,format_name(met["tag"], var_met, "", "pt" ,False))
                        met_phi = getattr(event,format_name(met["tag"], var_met, "", "phi",False))

                        # MET has to be corrected for different mu.pt definitions
                        if key_mu not in ["GEN", "PF"]:
                            met_px, met_py = met_pt*math.cos(met_phi),  met_pt*math.sin(met_phi)
                            for imuon,muon in enumerate(muons):
                                mu_pt = getattr(muon, format_name("",var_mu, "", "pt", False))
                                muon_phi = muon.phi
                                met_px -= (mu_pt*math.cos(muon_phi) - getattr(muon, format_name("","","","pt",False))*math.cos(muon_phi) )
                                met_py -= (mu_pt*math.sin(muon_phi) - getattr(muon, format_name("","","","pt",False))*math.sin(muon_phi))
                            #print "%s-%s before: %s, %s " % (key_mu, var_mu, met_pt,met_phi)
                            met_pt = math.sqrt(met_px**2 + met_py**2)
                            muon_phi = math.atan2(met_py,met_px)
                            #print "%s-%s after: %s, %s " % (key_mu, var_mu, met_pt,met_phi)

                        for imuon,muon in enumerate(muons) :
                            # Reco muons
                            if key_mu!="GEN":
                                mu_pt = getattr(muon, format_name("",var_mu, "", "pt", False))
                                muon_phi, muon_mass = muon.phi,muon.mass
                            # Gen muons
                            else:
                                genIdx = muon.genPartIdx
                                if genIdx>=0 and genIdx<len(genparticles):
                                    genMu = genparticles[genIdx]
                                    mu_pt, muon_phi, muon_mass = genMu.pt, genMu.phi, genMu.mass
                                else:
                                    mu_pt, muon_phi, muon_mass = 0.,0.,0.
                            (mt,hpt,hphi,met_par,met_per,h_par,h_per) = evaluateMt( mu_pt, muon_phi, muon_mass, met_pt, met_phi)
                            # mu_pfIso04 = getattr(muon, format_name("",var_mu, "", "pfRelIso04_all", False))*mu_pt

                            # W-like (filled only when there are 2 muons)
                            if event.Vtype in [2,3] and imuon in [event.Idx_mu1, event.Idx_mu2]:
                                metWlike_px, metWlike_py = met_pt*math.cos(met_phi), met_pt*math.sin(met_phi)
                                # nu-like muon is the OTHER muon wrt the one being analyzed
                                nuLike_muon = muons[event.Idx_mu2] if imuon==event.Idx_mu1 else muons[event.Idx_mu1]
                                if key_mu!="GEN":
                                    nuLike_mu_pt = getattr(nuLike_muon, format_name("",var_mu, "", "pt", False))
                                    nuLike_muon_phi,nuLike_muon_mass = nuLike_muon.phi,nuLike_muon.mass
                                else:
                                    nuLike_genIdx = nuLike_muon.genPartIdx
                                    if nuLike_genIdx>=0 and nuLike_genIdx<len(genparticles):
                                        nuLike_genMu = genparticles[nuLike_genIdx]
                                        nuLike_mu_pt, nuLike_muon_phi, nuLike_muon_mass = nuLike_genMu.pt, nuLike_genMu.phi, nuLike_genMu.mass
                                    else:
                                        nuLike_mu_pt, nuLike_muon_phi, nuLike_muon_mass = 0.,0.,0.
                                metWlike_px += nuLike_muon.pt*math.cos(nuLike_muon.phi)
                                metWlike_py += nuLike_muon.pt*math.sin(nuLike_muon.phi)
                                metWlike_pt = math.sqrt(metWlike_px**2 + metWlike_py**2)
                                metWlike_phi = math.atan2(metWlike_py, metWlike_px)
                                (Wlikemt,Wlikehpt,Wlikehphi,Wlikemet_par,Wlikemet_per,Wlikeh_par,Wlikeh_per) = evaluateMt( mu_pt, muon_phi, muon_mass, metWlike_pt , metWlike_phi)
                                Mt_Wlike_vec[imuon] = Wlikemt
                                Recoil_pt_Wlike_vec[imuon] = Wlikehpt
                                Recoil_phi_Wlike_vec[imuon] = Wlikehphi
                                met_par_Wlike_vec[imuon] = Wlikemet_par
                                met_per_Wlike_vec[imuon] = Wlikemet_per
                                Recoil_par_Wlike_vec[imuon] = Wlikeh_par
                                Recoil_per_Wlike_vec[imuon] = Wlikeh_per
                            Mt_vec[imuon] = mt
                            Recoil_pt_vec[imuon] = hpt
                            Recoil_phi_vec[imuon] = hphi
                            met_par_vec[imuon] = met_par
                            met_per_vec[imuon] = met_per
                            met_pt_vec[imuon] = met_pt
                            Recoil_par_vec[imuon] = h_par
                            Recoil_per_vec[imuon] = h_per


                        self.out.fillBranch("%s_mt" % branch_label_mu, Mt_vec)
                        if self.lightOutput<2:
                            self.out.fillBranch("%s_hpt"  % branch_label_mu, Recoil_pt_vec)
                            self.out.fillBranch("%s_hphi" % branch_label_mu, Recoil_phi_vec)
                            self.out.fillBranch("%s_uPar" % branch_label_mu, met_par_vec)
                            self.out.fillBranch("%s_uPer" % branch_label_mu, met_per_vec)
                            self.out.fillBranch("%s_metPt" % branch_label_mu, met_pt_vec)
                            self.out.fillBranch("%s_hPar" % branch_label_mu, Recoil_par_vec)
                            self.out.fillBranch("%s_hPer" % branch_label_mu, Recoil_per_vec)
                            # W-like variables
                            self.out.fillBranch("%s_Wlikemt" % branch_label_mu, Mt_Wlike_vec)
                            self.out.fillBranch("%s_Wlikehpt" % branch_label_mu, Recoil_pt_Wlike_vec)
                            self.out.fillBranch("%s_Wlikehphi" % branch_label_mu, Recoil_phi_Wlike_vec)
                            self.out.fillBranch("%s_WlikeuPar" % branch_label_mu, met_par_Wlike_vec)
                            self.out.fillBranch("%s_WlikeuPer" % branch_label_mu, met_per_Wlike_vec)
                            self.out.fillBranch("%s_WlikehPar" % branch_label_mu, Recoil_par_Wlike_vec)
                            self.out.fillBranch("%s_WlikehPer" % branch_label_mu, Recoil_per_Wlike_vec)

                        # MET projections along reco/gen Z
                        branch_label_Z  = ("GenFromRecoZ" if key_mu=="GEN" else "RecoZ")
                        branch_label_Z += format_name(mu["tag"], var_mu, met["tag"], var_met, True)
                        branch_label_Zmu  = ("GenFromRecoZ" if key_mu=="GEN" else "RecoZ")
                        branch_label_Zmu += format_name(mu["tag"], var_mu, "", "", True)
                        if event.Vtype in [2,3]:
                            Zpt = getattr(event, "%s_pt"  %  branch_label_Zmu)
                            Zphi = getattr(event, "%s_phi"  %  branch_label_Zmu)
                            Zmass = getattr(event, "%s_mass"  %  branch_label_Zmu)
                            (Zmt, Zh_pt, Zh_phi, Zmet_par, Zmet_per, Zh_par,Zh_per) = evaluateMt( Zpt, Zphi, Zmass, met_pt, met_phi)
                        else:
                            (Zmt, Zh_pt, Zh_phi, Zmet_par, Zmet_per,Zh_par,Zh_per) = (0., 0., 0., 0., 0., 0., 0.)
                        if self.lightOutput<2:
                            self.out.fillBranch("%s_uPar" % branch_label_Z, Zmet_par)
                            self.out.fillBranch("%s_uPer" % branch_label_Z, Zmet_per)
                            self.out.fillBranch("%s_hPar" % branch_label_Z, Zh_par)
                            self.out.fillBranch("%s_hPer" % branch_label_Z, Zh_per)

        '''
        for imuon,muon in enumerate(muons) : #definition of absolute isolation
            mu_pt = getattr(muon, "pt")
            mu_pfIso04 = getattr(muon, "pfRelIso04_all")*mu_pt
            pfIso04_vec[imuon] = mu_pfIso04
        self.out.fillBranch("Muon_pfIso04", pfIso04_vec)
        '''

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

addVarModule = lambda : additionalVariables()
