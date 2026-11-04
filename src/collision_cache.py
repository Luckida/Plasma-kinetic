"""Optional sparse cache of EXACTLY the same binary projection quadrature.

Only velocity-grid geometry is reused. Rates, screening and both distributions
are recomputed every call. Large grids fall back to the direct kernel within
the configured memory budget. The cache never contains atomic cross sections.
"""
import hashlib
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix,csr_matrix
import kinetic_kernels as K


class CollisionMaps:
    def __init__(self,pair,ge,gi,scatter_mu,scatter_w,scatter_phi,budget_mb,workers,root,progress=None):
        self.maps=[];self.workers=workers;self.pair=pair;self.mass_i=gi.mass
        ne,ni,naz=pair.rel.shape;raw_bytes=ne*ni*len(scatter_mu)*len(scatter_phi)*4*24
        budget=int(budget_mb*1024**2)
        if budget<=0 or raw_bytes>budget:return
        digest=hashlib.sha256(Path(K.__file__).read_bytes())
        for arr in (ge.e,ge.mu,gi.e,gi.mu,scatter_mu,scatter_w,scatter_phi):digest.update(np.asarray(arr,dtype=np.float64).tobytes())
        digest.update(str((naz,gi.mass)).encode());key=digest.hexdigest()[:24]
        folder=Path(root)/'.kinetic_cache';folder.mkdir(exist_ok=True)
        used=0
        for p in range(naz):
            if progress is not None:progress.module('Геометрия e–ion и компиляция Numba',f'азимут {p+1}/{naz}',force=True)
            if used+raw_bytes>budget:self.maps=[];return
            file=folder/f'{key}_{p}.npz'
            if file.exists():
                with np.load(file,allow_pickle=False) as z:
                    em=csr_matrix((z['ed'],z['ei'],z['ep']),shape=(ne*ne,ne*ni))
                    im=csr_matrix((z['id'],z['ii'],z['ip']),shape=(ne*ni,ne*ni))
            else:
                with K.kernel_threads(workers):
                    er,ew,ir,iw=K.scatter_map(pair.rel,pair.cx,pair.cy,pair.cz,p,ge.e,ge.mu,gi.e,gi.mu,scatter_mu,scatter_w,scatter_phi,gi.mass)
                columns=np.repeat(np.arange(ne*ni,dtype=np.int32),ew.shape[1])
                mask=ew.ravel()>0
                em=coo_matrix((ew.ravel()[mask],(er.ravel()[mask],columns[mask])),shape=(ne*ne,ne*ni)).tocsr()
                mask=iw.ravel()>0
                im=coo_matrix((iw.ravel()[mask],(ir.ravel()[mask],columns[mask])),shape=(ne*ni,ne*ni)).tocsr()
                del er,ew,ir,iw,columns,mask
                temp=file.with_suffix('.tmp.npz')
                np.savez_compressed(temp,ed=em.data,ei=em.indices,ep=em.indptr,id=im.data,ii=im.indices,ip=im.indptr)
                temp.replace(file)
            used+=sum(a.nbytes for a in (em.data,em.indices,em.indptr,im.data,im.indices,im.indptr))
            self.maps.append((em,im))
        self.memory_bytes=used

    def generator(self,ions,screen,logcap,ion_charge=1.):
        if not self.maps:raise ValueError('Sparse collision map unavailable')
        with K.kernel_threads(self.workers):
            rates=K.pair_rates(self.pair.rel,screen,logcap,self.mass_i,ion_charge)
        return self.generator_rates(ions,rates)

    def generator_rates(self,ions,rates):
        if not self.maps:raise ValueError('Sparse collision map unavailable')
        ne,ni,naz=self.pair.rel.shape;eg=np.zeros((ne,ne));ig=np.zeros((ne,ni));loss=np.zeros((ne,ni))
        with K.kernel_threads(self.workers):
            for p,(em,im) in enumerate(self.maps):
                weighted=(rates[:,:,p]*ions[None,:]/naz).ravel()
                eg+=K.sparse_product(em.data,em.indices,em.indptr,weighted).reshape(ne,ne)
                ig+=K.sparse_product(im.data,im.indices,im.indptr,weighted).reshape(ne,ni)
                loss+=weighted.reshape(ne,ni)
        return eg,ig,loss,loss.sum(axis=1)
