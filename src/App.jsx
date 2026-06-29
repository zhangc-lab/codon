import { useState, useCallback, useMemo, useEffect } from "react";

// ── Genetic code ──────────────────────────────────────────────────────────────
const GC = {TTT:'F',TTC:'F',TTA:'L',TTG:'L',CTT:'L',CTC:'L',CTA:'L',CTG:'L',ATT:'I',ATC:'I',ATA:'I',ATG:'M',GTT:'V',GTC:'V',GTA:'V',GTG:'V',TCT:'S',TCC:'S',TCA:'S',TCG:'S',AGT:'S',AGC:'S',CCT:'P',CCC:'P',CCA:'P',CCG:'P',ACT:'T',ACC:'T',ACA:'T',ACG:'T',GCT:'A',GCC:'A',GCA:'A',GCG:'A',TAT:'Y',TAC:'Y',TAA:'*',TAG:'*',TGA:'*',CAT:'H',CAC:'H',CAA:'Q',CAG:'Q',AAT:'N',AAC:'N',AAA:'K',AAG:'K',GAT:'D',GAC:'D',GAA:'E',GAG:'E',TGT:'C',TGC:'C',TGG:'W',CGT:'R',CGC:'R',CGA:'R',CGG:'R',AGA:'R',AGG:'R',GGT:'G',GGC:'G',GGA:'G',GGG:'G'};

// ── E. coli K12 (Kazusa) ──────────────────────────────────────────────────────
const EU_EC = {TTT:22.0,TTC:16.4,TTA:14.0,TTG:13.4,CTT:11.5,CTC:10.5,CTA:3.8,CTG:52.0,ATT:31.7,ATC:23.3,ATA:7.2,ATG:27.4,GTT:19.1,GTC:14.9,GTA:11.0,GTG:25.8,TCT:10.4,TCC:9.5,TCA:8.3,TCG:9.0,AGT:9.3,AGC:16.2,CCT:8.5,CCC:5.7,CCA:9.1,CCG:25.3,ACT:11.4,ACC:22.7,ACA:8.8,ACG:13.6,GCT:16.3,GCC:25.5,GCA:20.5,GCG:33.3,TAT:17.6,TAC:14.3,CAT:14.5,CAC:11.3,CAA:15.3,CAG:31.1,AAT:20.3,AAC:22.0,AAA:36.6,AAG:12.3,GAT:35.7,GAC:20.1,GAA:42.3,GAG:18.2,TGT:5.5,TGC:6.7,TGG:15.3,CGT:21.7,CGC:22.4,CGA:3.6,CGG:5.8,AGA:2.4,AGG:1.5,GGT:24.2,GGC:29.5,GGA:8.3,GGG:11.0};
const BC_EC = {F:'TTT',L:'CTG',I:'ATT',M:'ATG',V:'GTG',S:'AGC',P:'CCG',T:'ACC',A:'GCG',Y:'TAT',H:'CAT',Q:'CAG',N:'AAC',K:'AAA',D:'GAT',E:'GAA',C:'TGC',W:'TGG',R:'CGT',G:'GGC','*':'TAA'};

// ── Y. lipolytica CLIB89 (Kazusa) ─────────────────────────────────────────────
const EU_YL = {TTT:17.2,TTC:21.8,TTA:4.3,TTG:22.6,CTT:11.9,CTC:14.2,CTA:7.9,CTG:30.6,ATT:18.5,ATC:23.6,ATA:6.2,ATG:22.3,GTT:13.7,GTC:15.5,GTA:8.1,GTG:22.4,TCT:8.9,TCC:12.8,TCA:8.7,TCG:8.9,AGT:6.0,AGC:14.4,CCT:8.1,CCC:13.2,CCA:11.3,CCG:7.5,ACT:9.7,ACC:20.3,ACA:9.2,ACG:8.8,GCT:14.8,GCC:24.3,GCA:10.3,GCG:8.9,TAT:10.5,TAC:17.3,CAT:7.5,CAC:13.0,CAA:13.6,CAG:21.8,AAT:12.5,AAC:22.5,AAA:20.8,AAG:27.5,GAT:22.5,GAC:25.8,GAA:27.2,GAG:30.6,TGT:4.2,TGC:7.8,TGG:12.3,CGT:5.8,CGC:10.3,CGA:5.3,CGG:5.7,AGA:18.4,AGG:14.8,GGT:15.0,GGC:22.3,GGA:12.4,GGG:7.6};
const BC_YL = {F:'TTC',L:'CTG',I:'ATC',M:'ATG',V:'GTG',S:'AGC',P:'CCC',T:'ACC',A:'GCC',Y:'TAC',H:'CAC',Q:'CAG',N:'AAC',K:'AAG',D:'GAC',E:'GAG',C:'TGC',W:'TGG',R:'AGA',G:'GGC','*':'TAA'};

// ── Organism config ───────────────────────────────────────────────────────────
const ORGS = {
  ecoli:{
    label:'E. coli K12', eu:EU_EC, bc:BC_EC,
    strainRec:(vr,cp)=>[
      vr>0?{strain:'Rosetta2(DE3)',reason:'Very rare codons detected — supplemental tRNA plasmid required',color:'#f59e0b'}
          :{strain:'BL21(DE3)',reason:'No critically rare codons; standard expression strain sufficient',color:'#22c55e'},
      cp>0?{strain:'SHuffle T7 Express',reason:`${cp} close Cys pair(s) — cytoplasmic DsbC may improve folding`,color:'#f87171'}:null,
    ].filter(Boolean),
    cysNote:'Add reducing agent (1–5 mM DTT or 10 mM 2-ME) to all buffers. Verify any additional cofactor requirements specific to your protein.',
    footer:'E. coli K12 codon usage (Kazusa) · CAI: >0.8 good · >0.6 moderate',
  },
  ylip:{
    label:'Y. lipolytica', eu:EU_YL, bc:BC_YL,
    strainRec:(vr,cp)=>[
      vr>0?{strain:'Codon optimization recommended',reason:'Very rare Y. lipolytica codons present — full gene synthesis advised',color:'#f59e0b'}
          :{strain:'Po1g / W29',reason:'No critically rare codons; standard Y. lipolytica strain suitable',color:'#22c55e'},
      cp>0?{strain:'Verify disulfide bonds',reason:`${cp} close Cys pair(s) — Y. lipolytica secretory pathway supports disulfides; validate folding`,color:'#f87171'}:null,
    ].filter(Boolean),
    cysNote:'Y. lipolytica secretory pathway supports disulfide bond formation. Verify any additional cofactor and media requirements specific to your protein.',
    footer:'Y. lipolytica CLIB89 codon usage (Kazusa) · CAI: >0.8 good · >0.6 moderate',
  },
};

// ── Static data ───────────────────────────────────────────────────────────────
const AA_INFO = [
  {o:'A',t:'Ala',f:'Alanine'},{o:'R',t:'Arg',f:'Arginine'},{o:'N',t:'Asn',f:'Asparagine'},
  {o:'D',t:'Asp',f:'Aspartic acid'},{o:'C',t:'Cys',f:'Cysteine'},{o:'E',t:'Glu',f:'Glutamic acid'},
  {o:'Q',t:'Gln',f:'Glutamine'},{o:'G',t:'Gly',f:'Glycine'},{o:'H',t:'His',f:'Histidine'},
  {o:'I',t:'Ile',f:'Isoleucine'},{o:'L',t:'Leu',f:'Leucine'},{o:'K',t:'Lys',f:'Lysine'},
  {o:'M',t:'Met',f:'Methionine'},{o:'F',t:'Phe',f:'Phenylalanine'},{o:'P',t:'Pro',f:'Proline'},
  {o:'S',t:'Ser',f:'Serine'},{o:'T',t:'Thr',f:'Threonine'},{o:'W',t:'Trp',f:'Tryptophan'},
  {o:'Y',t:'Tyr',f:'Tyrosine'},{o:'V',t:'Val',f:'Valine'},
];
const CAI_COLOR = {Good:'#22c55e',Moderate:'#f59e0b',Poor:'#ef4444'};
const AA_COLORS = {C:'#f87171',M:'#fb923c',W:'#a78bfa',H:'#34d399',P:'#60a5fa',K:'#facc15',R:'#f472b6'};
const TABS = [{id:'rare',label:'Rare Codons'},{id:'clusters',label:'Stall Clusters'},{id:'cys',label:'Cys Analysis'},{id:'map',label:'Codon Map'},{id:'aa',label:'AA Composition'},{id:'opt',label:'Optimized Seq'}];

const EXAMPLES = [
  {label:'AgGPPS2',desc:'Terpene synthase',seq:'>AgGPPS2 (Abies grandis)\natgttcgatttcaacaaatacatggacagcaaagccatgaccgttaacgaagcgctgaacaaagcgatcccgcttcgttacccgcaaaaaatctacgaatccatgcgttactctctgctggcgggtgggaaacgcgttcgtccggtgctgtgtattgcggcgtgcgagctggtaggtggtaccgaagagctggctatcccgactgcatgcgcaatcgaaatgatccacaccatgtctctgatgcacgacgatctgccgtgtatcgacaacgacgacttacgtcgtggtaagccaactaaccacaagatctttggtgaagacacggctgtgaccgccggtaatgcactgcactcctacgcattcgaacacatcgcggtatctaccagcaaaactgtcggtgcggaccgtatcctgcgtatggtttctgaactgggtcgtgcaaccggtagcgaaggcgttatgggcggccagatggtagacattgcgtccgaaggcgatccgtccatcgatctgcagactctggaatggattcacattcataagaccgctatgttgctggagtgttccgtagtttgcggtgcaatcattggcggtgcgtctgaaattgtgatcgaacgtgctcgtcgttacgcgcgttgcgtggggctgctgttccaagtggtcgatgatatcttagatgttaccaaaagcagtgacgagctgggcaaaaccgcgggtaaagatctgatctctgacaaagcgacttacccaaaacttatgggcctggaaaaagcgaaagaattcagcgacgaactgctgaaccgtgctaagggcgaactgtcatgcttcgatccggtaaaagcagcgccactcctgggtctggctgactacgttgcattccgtcaaaactga'},
  {label:'eGFP',desc:'Fluorescent reporter',seq:'>eGFP\natgcatcatcatcaccatcacgagctcATGgtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagtaa'},
  {label:'ispA',desc:'FPP synthase (E. coli)',seq:'>ispA\natgcatcatcatcaccatcacgagctcgactttccgcagcaactcgaagcctgcgttaagcaggccaaccaggcgctgagccgttttatcgccccactgccctttcagaacactcccgtggtcgaaaccatgcagtatggcgcattattaggtggtaagcgcctgcgacctttcctggtttatgccaccggtcatatgttcggcgttagcacaaacacgctggacgcacccgctgccgccgttgagtgtatccacgcttactcattaattcatgatgatttaccggcaatggatgatgacgatctgcgtcgcggtttgccaacctgccatgtgaagtttggcgaagcaaacgcgattctcgctggcgacgctttacaaacgctggcgttctcgattttaagcgatgccgatatgccggaagtgtcggaccgcgacagaatttcgatgatttctgaactggcgagcgccagtggtattgccggaatgtgcggtggtcaggcattagatttagacgcggaaggcaaacacgtacctctggacgcgcttgagcgtattcatcgtcataaaaccggcgcattgattcgcgccgccgttcgccttggtgcattaagcgccggagataaaggacgtcgtgctctgccggtactcgacaagtatgcagagagcatcggccttgccttccaggttcaggatgacatcctggatgtggtgggagatactgcaacgttgggaaaacgccagggtgccgaccagcaacttggtaaaagtacctaccctgcacttctgggtcttgagcaagcccggaagaaagcccgggatctgatcgacgatgcccgtcagtcgctgaaacaactggctgaacagtcactcgatacctcggcactggaagcgctagcggactacatcatccagcgtaataaataa'},
];

// ── FASTA parser ──────────────────────────────────────────────────────────────
function parseFASTA(raw) {
  const trimmed=raw.trim();
  if(!trimmed) return [];
  if(!trimmed.includes('>')) {
    const seq=trimmed.toUpperCase().replace(/[^ATGCN\s]/g,'').replace(/\s/g,'');
    return seq?[{name:'Sequence 1',seq}]:[];
  }
  const seqs=[]; let cur=null;
  for(const line of trimmed.split('\n')) {
    const t=line.trim();
    if(t.startsWith('>')) {if(cur?.seq)seqs.push(cur);cur={name:t.slice(1).trim()||'Unnamed',seq:''};}
    else if(t&&cur) cur.seq+=t.toUpperCase().replace(/[^ATGCN]/g,'');
  }
  if(cur?.seq) seqs.push(cur);
  return seqs;
}

// ── Analysis engine ───────────────────────────────────────────────────────────
function analyze(rawSeq, eu, bc) {
  let seq=rawSeq.toUpperCase().replace(/[^ATGCN]/g,'');
  if(seq.length%3!==0) seq=seq.slice(0,seq.length-(seq.length%3));
  const codons=[];
  for(let i=0;i<seq.length-2;i+=3) codons.push(seq.slice(i,i+3));
  const protein=codons.map(c=>GC[c]||'?').join('');
  const aas=protein.replace(/\*/g,'').split('');
  const rareCodons=codons.map((c,i)=>{
    const aa=GC[c];if(!aa||aa==='*') return null;
    const freq=eu[c]??0;
    return freq<10?{position:i+1,codon:c,aa,freq,veryRare:freq<5}:null;
  }).filter(Boolean);
  const mask=codons.map(c=>(eu[c]??99)<10);
  const clusters=[]; let i=0;
  while(i<mask.length){
    if(mask[i]){let j=i;while(j<mask.length&&mask[j])j++;
      if(j-i>=2)clusters.push({start:i+1,end:j,codons:codons.slice(i,j).map((c,k)=>({pos:i+k+1,codon:c,aa:GC[c]||'?',freq:eu[c]??0}))});i=j;}else i++;}
  const replacements=rareCodons.map(r=>{
    const opt=bc[r.aa]||r.codon;if(opt===r.codon) return null;
    return {...r,optimized:opt,optFreq:eu[opt]??0};
  }).filter(Boolean);
  const optCodons=codons.map(c=>{const aa=GC[c];return aa&&aa!=='*'&&(eu[c]??0)<10?(bc[aa]||c):c;});
  const optimizedSeq=optCodons.join('');
  const caiVals=codons.filter(c=>GC[c]&&GC[c]!=='*').map(c=>{
    const aa=GC[c];return Math.log(Math.max(eu[c]||0.01,0.01)/Math.max(eu[bc[aa]||c]||1,0.01));});
  const cai=caiVals.length?Math.exp(caiVals.reduce((a,b)=>a+b,0)/caiVals.length):0;
  const cysPos=codons.map((c,i)=>GC[c]==='C'?i+1:null).filter(Boolean);
  const closePairs=[];
  for(let a=0;a<cysPos.length;a++)for(let b=a+1;b<cysPos.length;b++)
    if(cysPos[b]-cysPos[a]<=6)closePairs.push({pos1:cysPos[a],pos2:cysPos[b],dist:cysPos[b]-cysPos[a]});
  const aaCounts={};
  for(const aa of aas) aaCounts[aa]=(aaCounts[aa]||0)+1;
  return {
    summary:{lengthBp:seq.length,proteinLength:aas.length,
      cai:Math.round(cai*1000)/1000,caiGrade:cai>0.8?'Good':cai>0.6?'Moderate':'Poor',
      rareCount:rareCodons.length,veryRareCount:rareCodons.filter(r=>r.veryRare).length,clusterCount:clusters.length},
    protein,originalSeq:seq,optimizedSeq,rareCodons,clusters,replacements,aaCounts,bc,
    cys:{count:cysPos.length,positions:cysPos,closePairs},
  };
}

// ── Export helpers ────────────────────────────────────────────────────────────
function dl(content,filename,type){
  const url=URL.createObjectURL(new Blob([content],{type}));
  Object.assign(document.createElement('a'),{href:url,download:filename}).click();
  URL.revokeObjectURL(url);
}
function toMD(name,r,orgLabel){
  const s=r.summary;
  return [
    `# Codon Analysis: ${name}`,`*${new Date().toLocaleString()} | ${orgLabel}*\n`,
    '## Summary','| Parameter | Value |','|-----------|-------|',
    `| Reference organism | ${orgLabel} |`,
    `| Length | ${s.lengthBp} bp |`,`| Protein | ${s.proteinLength} aa |`,
    `| CAI | ${s.cai} (${s.caiGrade}) |`,`| Rare Codons | ${s.rareCount} (${s.veryRareCount} very rare) |`,
    `| Stall Clusters | ${s.clusterCount} |`,`| Cys | ${r.cys.count} residues, ${r.cys.closePairs.length} close pair(s) |\n`,
    '## Rare Codons',
    r.replacements.length
      ?['| Pos | Codon | AA | Freq | → Opt | New Freq |','|-----|-------|-----|-----:|-------|--------:|',
        ...r.replacements.map(x=>`| ${x.position} | ${x.codon}${x.veryRare?' ⚠':''} | ${x.aa} | ${x.freq} | ${x.optimized} | ${x.optFreq} |`)].join('\n')
      :'No rare codons found.',
    '\n## Stall Clusters',
    r.clusters.length?r.clusters.map(cl=>`**Pos ${cl.start}–${cl.end}:** ${cl.codons.map(c=>`${c.codon}(${c.aa})`).join(' ')}`).join('\n'):'None.',
    '\n## Cysteine','Positions: '+(r.cys.positions.join(', ')||'None'),
    r.cys.closePairs.map(p=>`- Cys${p.pos1}–Cys${p.pos2} (${p.dist} aa)`).join('\n')||'No close pairs.',
    '\n## Amino Acid Composition','| Amino Acid | 3-Letter | 1-Letter | Count | % |','|------------|----------|----------|------:|--:|',
    ...AA_INFO.map(a=>{const n=r.aaCounts[a.o]||0;return `| ${a.f} | ${a.t} | ${a.o} | ${n} | ${s.proteinLength?(n/s.proteinLength*100).toFixed(1):0}% |`;}),
    '\n## Optimized Sequence','```\n'+r.optimizedSeq+'\n```',
    '\n## Protein Sequence','```\n'+r.protein.replace(/\*$/,'')+'\n```',
  ].join('\n');
}
function toWordHTML(name,r,orgLabel){
  const s=r.summary,cc=CAI_COLOR[s.caiGrade]||'#666';
  const aaMax=Math.max(...AA_INFO.map(a=>r.aaCounts[a.o]||0),1);
  return `<!DOCTYPE html><html xmlns:o="urn:schemas-microsoft-com:office:office" xmlns:w="urn:schemas-microsoft-com:office:word">
<head><meta charset="UTF-8"><style>
body{font-family:Calibri,Arial,sans-serif;font-size:11pt;margin:2cm}
h1{color:#1e3a5f;border-bottom:2px solid #1e3a5f;padding-bottom:6px}
h2{color:#1e3a5f;border-bottom:1px solid #ccc;margin-top:20px}
table{border-collapse:collapse;width:100%;margin:10px 0}
th{background:#1e3a5f;color:#fff;padding:6px 9px;text-align:left}
td{border:1px solid #ddd;padding:5px 9px}tr:nth-child(even){background:#f8fafc}
pre{background:#f5f5f5;padding:10px;font-family:'Courier New',monospace;font-size:9pt;word-break:break-all;white-space:pre-wrap;border:1px solid #ddd}
</style></head><body>
<h1>Codon Analysis: ${name}</h1><p><em>${new Date().toLocaleString()} | ${orgLabel}</em></p>
<h2>Summary</h2><table><thead><tr><th>Parameter</th><th>Value</th></tr></thead><tbody>
<tr><td>Reference organism</td><td>${orgLabel}</td></tr>
<tr><td>Length</td><td>${s.lengthBp} bp</td></tr><tr><td>Protein</td><td>${s.proteinLength} aa</td></tr>
<tr><td>CAI</td><td style="color:${cc};font-weight:bold">${s.cai} (${s.caiGrade})</td></tr>
<tr><td>Rare Codons</td><td>${s.rareCount} (${s.veryRareCount} very rare)</td></tr>
<tr><td>Stall Clusters</td><td>${s.clusterCount}</td></tr>
<tr><td>Cys Residues</td><td>${r.cys.count} (${r.cys.closePairs.length} close pair(s))</td></tr>
</tbody></table>
<h2>Rare Codons</h2>
${r.replacements.length
  ?`<table><thead><tr><th>Pos</th><th>Codon</th><th>AA</th><th>Freq/1000</th><th>→ Opt</th><th>New Freq</th></tr></thead><tbody>${r.replacements.map(x=>`<tr><td>${x.position}</td><td style="color:${x.veryRare?'#dc2626':'#d97706'};font-weight:bold">${x.codon}${x.veryRare?' ⚠':''}</td><td>${x.aa}</td><td>${x.freq}</td><td style="color:#16a34a;font-weight:bold">${x.optimized}</td><td>${x.optFreq}</td></tr>`).join('')}</tbody></table>`
  :'<p style="color:#16a34a">✓ No rare codons found.</p>'}
<h2>Stall Clusters</h2>
${r.clusters.length?r.clusters.map(cl=>`<div style="background:#fff5f5;border:1px solid #fca5a5;padding:8px;margin:6px 0"><b style="color:#dc2626">⚠ Pos ${cl.start}–${cl.end}:</b> ${cl.codons.map(c=>`<code>${c.codon}(${c.aa},${c.freq})</code>`).join(' ')}</div>`).join(''):'<p style="color:#16a34a">✓ No clusters.</p>'}
<h2>Cysteine Analysis</h2><p>Positions: ${r.cys.positions.join(', ')||'None'}</p>
${r.cys.closePairs.length?`<p style="color:#dc2626">⚠ Close pairs (≤6 aa):</p><ul>${r.cys.closePairs.map(p=>`<li>Cys${p.pos1}–Cys${p.pos2} (${p.dist} aa apart)</li>`).join('')}</ul>`:'<p style="color:#16a34a">✓ No close Cys pairs.</p>'}
<h2>Amino Acid Composition</h2>
<table><thead><tr><th>Amino Acid</th><th>3-Letter</th><th>1-Letter</th><th>Count</th><th>Bar</th><th>%</th></tr></thead><tbody>
${AA_INFO.map(a=>{const n=r.aaCounts[a.o]||0,pct=s.proteinLength?(n/s.proteinLength*100).toFixed(1):0;
return `<tr ${n===0?'style="color:#aaa"':''}><td>${a.f}</td><td>${a.t}</td><td style="font-family:monospace;font-weight:bold">${a.o}</td><td>${n}</td><td><div style="background:#ddd;height:8px;border-radius:3px"><div style="background:#1d4ed8;height:8px;border-radius:3px;width:${Math.round(n/aaMax*100)}%"></div></div></td><td>${pct}%</td></tr>`;}).join('')}
</tbody></table>
<h2>Optimized Sequence (${r.replacements.length} codon(s) changed)</h2><pre>${r.optimizedSeq}</pre>
<h2>Protein Sequence</h2><pre>${r.protein.replace(/\*$/,'')}</pre>
</body></html>`;
}

// ── Sub-components ────────────────────────────────────────────────────────────
function CAIGauge({cai,grade}){
  const pct=Math.round(cai*100),color=CAI_COLOR[grade]||'#94a3b8';
  return(
    <div style={{display:'flex',flexDirection:'column',alignItems:'center',gap:4}}>
      <svg width={120} height={70} viewBox="0 0 120 70">
        <path d="M10 60 A50 50 0 0 1 110 60" fill="none" stroke="#1e293b" strokeWidth={12} strokeLinecap="round"/>
        <path d="M10 60 A50 50 0 0 1 110 60" fill="none" stroke={color} strokeWidth={12} strokeLinecap="round" strokeDasharray={`${pct*1.571} 157.1`}/>
        <text x={60} y={58} textAnchor="middle" fill={color} fontSize={20} fontWeight={700}>{cai}</text>
      </svg>
      <span style={{fontSize:12,color,fontWeight:600,textTransform:'uppercase'}}>{grade}</span>
    </div>
  );
}

function CodonMap({seq,replacements}){
  if(!seq) return null;
  const rs={};replacements.forEach(r=>{rs[r.position-1]=r;});
  const codons=[];for(let i=0;i<seq.length;i+=3) codons.push(seq.slice(i,i+3));
  return(
    <div style={{fontFamily:'monospace',fontSize:11,lineHeight:'2.1',background:'#0a0f1a',borderRadius:8,padding:'12px 14px',border:'1px solid #1e293b',display:'flex',flexWrap:'wrap',gap:'2px 4px',maxHeight:320,overflowY:'auto'}}>
      {codons.map((c,i)=>{const r=rs[i];
        return <span key={i} title={r?`${r.aa} — ${r.freq}/1000 → ${r.optimized} (${r.optFreq}/1000)`:undefined}
          style={{background:r?.veryRare?'#7f1d1d':r?'#78350f':'transparent',color:r?.veryRare?'#fca5a5':r?'#fde68a':'#475569',borderRadius:2,padding:'0 2px',borderBottom:r?`2px solid ${r.veryRare?'#ef4444':'#f59e0b'}`:'2px solid transparent',cursor:r?'help':'default',display:'inline-block'}}>{c}</span>;
      })}
    </div>
  );
}

function OptimizedDisplay({originalSeq,optimizedSeq}){
  const orig=[],opt=[];
  for(let i=0;i<originalSeq.length;i+=3){orig.push(originalSeq.slice(i,i+3));opt.push(optimizedSeq.slice(i,i+3));}
  const changed=orig.filter((c,i)=>c!==opt[i]).length;
  return(
    <div>
      <div style={{fontSize:11,color:'#64748b',marginBottom:8,display:'flex',gap:16,alignItems:'center',flexWrap:'wrap'}}>
        <span style={{color:'#94a3b8'}}><strong style={{color:'#f1f5f9'}}>{changed}</strong> codon(s) changed</span>
        <span style={{display:'flex',alignItems:'center',gap:5}}>
          <span style={{background:'#451a03',color:'#fbbf24',padding:'0 5px',borderRadius:2,fontFamily:'monospace',fontWeight:700}}>A</span>Changed nucleotide
        </span>
        <span style={{display:'flex',alignItems:'center',gap:5}}>
          <span style={{borderBottom:'2px solid #f59e0b',padding:'0 5px',fontFamily:'monospace',color:'#475569'}}>ATG</span>Changed codon
        </span>
      </div>
      <div style={{fontFamily:'monospace',fontSize:11,lineHeight:'2.1',background:'#0a0f1a',borderRadius:8,padding:'12px 14px',border:'1px solid #1e293b',display:'flex',flexWrap:'wrap',gap:'2px 4px',maxHeight:320,overflowY:'auto'}}>
        {opt.map((oc,ci)=>{
          const ow=orig[ci]||'',cc=oc!==ow;
          return(
            <span key={ci} style={{display:'inline-block',borderRadius:2,padding:'0 1px',borderBottom:cc?'2px solid #f59e0b':'2px solid transparent'}}>
              {oc.split('').map((nt,ni)=>{const chg=cc&&nt!==ow[ni];return <span key={ni} style={{color:chg?'#fbbf24':'#475569',background:chg?'#451a03':'transparent',borderRadius:1}}>{nt}</span>;})}
            </span>
          );
        })}
      </div>
    </div>
  );
}

function CopyBtn({text,label='Copy'}){
  const [done,setDone]=useState(false);
  return <button onClick={()=>{navigator.clipboard.writeText(text);setDone(true);setTimeout(()=>setDone(false),2000);}}
    style={{padding:'4px 10px',fontSize:11,background:done?'#166534':'#1e293b',color:done?'#86efac':'#94a3b8',border:'1px solid #334155',borderRadius:4,cursor:'pointer'}}>
    {done?'✓ Copied':label}</button>;
}

function SortArrow({dir}){
  return <span style={{fontSize:9,marginLeft:3,color:dir?'#38bdf8':'#e2e8f0'}}>{dir==='asc'?'▲':'▼'}</span>;
}

// ── Main App ──────────────────────────────────────────────────────────────────
export default function CodonAnalyzer(){
  const [input,setInput]=useState('');
  const [sequences,setSequences]=useState([]);
  const [selected,setSelected]=useState(0);
  const [error,setError]=useState('');
  const [tab,setTab]=useState('rare');
  const [sortAA,setSortAA]=useState('name-asc');
  const [org,setOrg]=useState('ecoli');

  const runAnalysis=useCallback((seqInput,orgKey)=>{
    const parsed=parseFASTA(seqInput);
    if(!parsed.length) return null;
    const {eu,bc}=ORGS[orgKey];
    return parsed.map(p=>({name:p.name,seq:p.seq,result:analyze(p.seq,eu,bc)}));
  },[]);

  const run=useCallback(()=>{
    setError('');
    const parsed=parseFASTA(input);
    if(!parsed.length){setError('No valid sequence found.');return;}
    for(const p of parsed) if(p.seq.length<9){setError(`"${p.name}" is too short (min 9 bp).`);return;}
    try{const res=runAnalysis(input,org);if(res){setSequences(res);setSelected(0);setTab('rare');}}
    catch(e){setError('Analysis failed: '+e.message);}
  },[input,org,runAnalysis]);

  // Re-analyze on organism switch if sequences are already loaded
  useEffect(()=>{
    if(sequences.length>0&&input.trim()){
      try{const res=runAnalysis(input,org);if(res)setSequences(res);}catch(e){}
    }
  },[org]); // eslint-disable-line

  const cur=sequences[selected],r=cur?.result,s=r?.summary;
  const cColor=s?CAI_COLOR[s.caiGrade]:'#94a3b8';
  const orgCfg=ORGS[org];

  const sortedAA=useMemo(()=>{
    if(!r) return AA_INFO;
    return [...AA_INFO].sort((a,b)=>{
      if(sortAA==='name-asc') return a.f.localeCompare(b.f);
      if(sortAA==='name-desc') return b.f.localeCompare(a.f);
      if(sortAA==='count-desc') return (r.aaCounts[b.o]||0)-(r.aaCounts[a.o]||0);
      if(sortAA==='count-asc') return (r.aaCounts[a.o]||0)-(r.aaCounts[b.o]||0);
      return 0;
    });
  },[r,sortAA]);

  const toggleSort=col=>{
    if(col==='name') setSortAA(s=>s==='name-asc'?'name-desc':'name-asc');
    if(col==='count') setSortAA(s=>s==='count-desc'?'count-asc':'count-desc');
  };
  const nameDir=sortAA.startsWith('name')?sortAA.split('-')[1]:null;
  const countDir=sortAA.startsWith('count')?sortAA.split('-')[1]:null;

  return(
    <div style={{minHeight:'100vh',background:'#0a0f1a',color:'#e2e8f0',fontFamily:"'Inter',system-ui,sans-serif",padding:'20px 16px'}}>
      {/* Header */}
      <div style={{maxWidth:920,margin:'0 auto 20px'}}>
        <div style={{display:'flex',alignItems:'center',gap:10,marginBottom:4}}>
          <svg width={28} height={28} viewBox="0 0 28 28" fill="none">
            <rect width={28} height={28} rx={6} fill="#1e3a5f"/>
            {[[8,9,'#38bdf8'],[14,14,'#22c55e'],[20,9,'#f59e0b'],[8,19,'#f59e0b'],[20,19,'#38bdf8']].map(([cx,cy,f],i)=><circle key={i} cx={cx} cy={cy} r={2.5} fill={f}/>)}
            {[[8,9,14,14],[20,9,14,14],[8,19,14,14],[20,19,14,14]].map(([x1,y1,x2,y2],i)=><line key={i} x1={x1} y1={y1} x2={x2} y2={y2} stroke="#334155" strokeWidth={1.5}/>)}
          </svg>
          <h1 style={{margin:0,fontSize:20,fontWeight:700,color:'#f1f5f9'}}>Codon Analyzer</h1>
          {/* Organism toggle */}
          <div style={{marginLeft:'auto',display:'flex',border:'1px solid #334155',borderRadius:6,overflow:'hidden'}}>
            {Object.entries(ORGS).map(([key,o])=>(
              <button key={key} onClick={()=>setOrg(key)} style={{
                padding:'5px 13px',fontSize:11,border:'none',cursor:'pointer',fontWeight:org===key?600:400,
                background:org===key?'#1d4ed8':'#1e293b',color:org===key?'#fff':'#64748b',whiteSpace:'nowrap',
                borderRight:key==='ecoli'?'1px solid #334155':'none',transition:'background 0.15s',
              }}>{o.label}</button>
            ))}
          </div>
        </div>
        <p style={{margin:0,fontSize:13,color:'#64748b'}}>Multi-gene · Rare codon detection · CAI · Cys risk · Optimization · Export</p>
      </div>

      <div style={{maxWidth:920,margin:'0 auto',display:'flex',flexDirection:'column',gap:16}}>
        {/* Input */}
        <div style={{background:'#0f172a',border:'1px solid #1e293b',borderRadius:10,padding:14}}>
          <div style={{display:'flex',justifyContent:'space-between',alignItems:'center',marginBottom:8}}>
            <span style={{fontSize:12,color:'#64748b',textTransform:'uppercase',letterSpacing:'0.05em'}}>DNA Sequence(s) — raw or multi-gene FASTA</span>
            {input&&<button onClick={()=>{setInput('');setSequences([]);}} style={{fontSize:11,color:'#475569',background:'transparent',border:'none',cursor:'pointer'}}>Clear</button>}
          </div>
          <textarea value={input} onChange={e=>setInput(e.target.value)}
            placeholder={">Gene1\nATGTTCGATTTCAACAAATACATGG...\n\n>Gene2\nATGGCGACCATTGCGATT..."}
            style={{width:'100%',minHeight:110,background:'#0a0f1a',color:'#94a3b8',border:'1px solid #1e293b',borderRadius:6,padding:'8px 10px',fontFamily:'monospace',fontSize:12,resize:'vertical',outline:'none',boxSizing:'border-box',lineHeight:1.6}}/>
          {error&&<p style={{margin:'6px 0 0',color:'#f87171',fontSize:12}}>⚠ {error}</p>}
          <div style={{marginTop:10,display:'flex',gap:8,flexWrap:'wrap'}}>
            <button onClick={run} style={{padding:'8px 22px',background:'#1d4ed8',color:'#fff',border:'none',borderRadius:6,fontWeight:600,fontSize:13,cursor:'pointer'}}>Analyze →</button>
            <button onClick={()=>setInput(EXAMPLES.map(e=>e.seq).join('\n\n'))}
              style={{padding:'6px 14px',background:'#172554',color:'#93c5fd',border:'1px solid #1d4ed8',borderRadius:6,fontSize:12,cursor:'pointer',fontWeight:600}}>
              Load the example
            </button>
          </div>
        </div>

        {/* Gene selector */}
        {sequences.length>1&&(
          <div style={{background:'#0f172a',border:'1px solid #1e293b',borderRadius:8,padding:'10px 14px'}}>
            <div style={{fontSize:11,color:'#475569',marginBottom:8,textTransform:'uppercase',letterSpacing:'0.05em'}}>{sequences.length} sequences — select to view:</div>
            <div style={{display:'flex',flexWrap:'wrap',gap:6}}>
              {sequences.map((seq,i)=>(
                <button key={i} onClick={()=>{setSelected(i);setTab('rare');}}
                  style={{padding:'5px 14px',fontSize:12,fontWeight:selected===i?600:400,background:selected===i?'#1d4ed8':'#1e293b',color:selected===i?'#fff':'#94a3b8',border:`1px solid ${selected===i?'#3b82f6':'#334155'}`,borderRadius:6,cursor:'pointer',display:'flex',alignItems:'center',gap:6}}>
                  {seq.name.length>26?seq.name.slice(0,24)+'…':seq.name}
                  {seq.result.summary.veryRareCount>0&&<span style={{background:'#7f1d1d',color:'#fca5a5',fontSize:9,padding:'1px 4px',borderRadius:2}}>{seq.result.summary.veryRareCount}⚠</span>}
                </button>
              ))}
            </div>
          </div>
        )}

        {r&&cur&&<>
          {sequences.length===1&&<div style={{fontSize:13,color:'#475569'}}>Analyzing: <strong style={{color:'#94a3b8'}}>{cur.name}</strong> <span style={{fontSize:11,color:'#334155',marginLeft:6}}>ref: {orgCfg.label}</span></div>}

          {/* Summary cards */}
          <div style={{display:'grid',gridTemplateColumns:'repeat(auto-fit,minmax(130px,1fr))',gap:10}}>
            {[
              {label:'Length',value:`${s.lengthBp} bp`,sub:`${s.proteinLength} aa`},
              {label:'CAI Score',value:s.cai,sub:s.caiGrade,color:cColor},
              {label:'Rare Codons',value:s.rareCount,sub:`${s.veryRareCount} very rare (<5)`,color:s.veryRareCount>0?'#f87171':s.rareCount>0?'#f59e0b':'#22c55e'},
              {label:'Stall Clusters',value:s.clusterCount,sub:'consecutive rare',color:s.clusterCount>0?'#f87171':'#22c55e'},
              {label:'Cys Residues',value:r.cys.count,sub:`${r.cys.closePairs.length} close pair(s)`,color:r.cys.closePairs.length>0?'#f87171':'#94a3b8'},
            ].map(card=>(
              <div key={card.label} style={{background:'#0f172a',border:'1px solid #1e293b',borderRadius:8,padding:'12px 14px'}}>
                <div style={{fontSize:11,color:'#475569',letterSpacing:'0.05em',textTransform:'uppercase',marginBottom:4}}>{card.label}</div>
                <div style={{fontSize:22,fontWeight:700,color:card.color||'#f1f5f9',lineHeight:1.1}}>{card.value}</div>
                <div style={{fontSize:11,color:'#475569',marginTop:2}}>{card.sub}</div>
              </div>
            ))}
          </div>

          {/* CAI + Strain */}
          <div style={{display:'flex',gap:12,flexWrap:'wrap'}}>
            <div style={{background:'#0f172a',border:'1px solid #1e293b',borderRadius:8,padding:'14px 18px',display:'flex',flexDirection:'column',alignItems:'center',gap:6,minWidth:160}}>
              <div style={{fontSize:11,color:'#475569',textTransform:'uppercase',letterSpacing:'0.05em'}}>CAI Gauge</div>
              <CAIGauge cai={s.cai} grade={s.caiGrade}/>
            </div>
            <div style={{flex:1,background:'#0f172a',border:'1px solid #1e293b',borderRadius:8,padding:'14px 16px',minWidth:240}}>
              <div style={{display:'flex',alignItems:'center',gap:8,marginBottom:10}}>
                <div style={{fontSize:11,color:'#475569',textTransform:'uppercase',letterSpacing:'0.05em'}}>Strain / Expression Recommendation</div>
                <span style={{fontSize:10,color:'#475569',background:'#1e293b',padding:'1px 7px',borderRadius:10,border:'1px solid #334155'}}>{orgCfg.label}</span>
              </div>
              {orgCfg.strainRec(s.veryRareCount,r.cys.closePairs.length).map((rec,i)=>(
                <div key={i} style={{display:'flex',gap:8,alignItems:'flex-start',marginBottom:8}}>
                  <span style={{color:rec.color,fontSize:14,marginTop:1}}>●</span>
                  <div><span style={{fontWeight:600,fontSize:13,color:rec.color}}>{rec.strain}</span>
                    <div style={{fontSize:12,color:'#64748b',lineHeight:1.4}}>{rec.reason}</div>
                  </div>
                </div>
              ))}
              {r.cys.closePairs.length>0&&<div style={{fontSize:11,color:'#78716c',background:'#1c1917',borderRadius:4,padding:'6px 8px',marginTop:4,borderLeft:'3px solid #78350f'}}>
                ⚠ {orgCfg.cysNote}
              </div>}
            </div>
          </div>

          {/* Tabs */}
          <div style={{background:'#0f172a',border:'1px solid #1e293b',borderRadius:10,overflow:'hidden'}}>
            <div style={{display:'flex',borderBottom:'1px solid #1e293b',overflowX:'auto'}}>
              {TABS.map(t=>(
                <button key={t.id} onClick={()=>setTab(t.id)} style={{padding:'10px 16px',fontSize:12,fontWeight:tab===t.id?600:400,color:tab===t.id?'#38bdf8':'#64748b',background:'transparent',border:'none',borderBottom:`2px solid ${tab===t.id?'#38bdf8':'transparent'}`,cursor:'pointer',whiteSpace:'nowrap'}}>{t.label}</button>
              ))}
            </div>
            <div style={{padding:16}}>

              {tab==='rare'&&(r.replacements.length===0
                ?<div style={{color:'#22c55e',fontSize:13}}>✓ No rare codons — sequence is well-adapted for <em>{orgCfg.label}</em>.</div>
                :<div style={{overflowX:'auto'}}>
                  <table style={{width:'100%',borderCollapse:'collapse',fontSize:12}}>
                    <thead><tr style={{color:'#475569',textAlign:'left'}}>
                      {['Position','Codon','AA','Freq/1000','→ Optimized','New Freq','Δ'].map(h=>(
                        <th key={h} style={{padding:'6px 10px',borderBottom:'1px solid #1e293b',fontWeight:600,fontSize:11,whiteSpace:'nowrap'}}>{h}</th>))}
                    </tr></thead>
                    <tbody>{r.replacements.map((x,i)=>(
                      <tr key={i} style={{borderBottom:'1px solid #1e293b',background:i%2?'transparent':'#0a0f1a'}}>
                        <td style={{padding:'6px 10px',color:'#64748b',fontFamily:'monospace'}}>{x.position}</td>
                        <td style={{padding:'6px 10px',fontFamily:'monospace',color:x.veryRare?'#f87171':'#fde68a',fontWeight:600}}>{x.codon}{x.veryRare&&<span style={{fontSize:10}}> ⚠</span>}</td>
                        <td style={{padding:'6px 10px',fontFamily:'monospace',color:'#94a3b8'}}>{x.aa}</td>
                        <td style={{padding:'6px 10px',color:x.veryRare?'#f87171':'#f59e0b'}}>{x.freq}</td>
                        <td style={{padding:'6px 10px',fontFamily:'monospace',color:'#22c55e',fontWeight:600}}>{x.optimized}</td>
                        <td style={{padding:'6px 10px',color:'#22c55e'}}>{x.optFreq}</td>
                        <td style={{padding:'6px 10px',color:'#22c55e',fontSize:11}}>+{(x.optFreq-x.freq).toFixed(1)}</td>
                      </tr>))}
                    </tbody>
                  </table>
                </div>
              )}

              {tab==='clusters'&&(r.clusters.length===0
                ?<div style={{color:'#22c55e',fontSize:13}}>✓ No consecutive rare codon clusters found.</div>
                :r.clusters.map((cl,i)=>(
                  <div key={i} style={{background:'#1c0a0a',border:'1px solid #7f1d1d',borderRadius:6,padding:'10px 12px',marginBottom:8}}>
                    <div style={{fontSize:12,fontWeight:600,color:'#fca5a5',marginBottom:6}}>⚠ Cluster: positions {cl.start}–{cl.end} <span style={{color:'#78716c',fontWeight:400}}>(ribosome stall risk)</span></div>
                    <div style={{display:'flex',gap:8,flexWrap:'wrap'}}>
                      {cl.codons.map((c,j)=>(<div key={j} style={{background:'#450a0a',borderRadius:4,padding:'4px 8px'}}>
                        <span style={{fontFamily:'monospace',color:'#fca5a5',fontWeight:600}}>{c.codon}</span>
                        <span style={{color:'#78716c',fontSize:10,marginLeft:4}}>({c.aa},{c.freq})</span>
                      </div>))}
                    </div>
                    <div style={{fontSize:11,color:'#64748b',marginTop:6}}>Replace: {cl.codons.map(c=>r.bc[c.aa]||c.codon).join(' + ')}</div>
                  </div>
                ))
              )}

              {tab==='cys'&&(
                <div>
                  <div style={{display:'flex',gap:20,flexWrap:'wrap',marginBottom:12}}>
                    <div style={{fontSize:13,color:'#94a3b8'}}><strong style={{color:'#f1f5f9'}}>{r.cys.count}</strong> Cys residues</div>
                    <div style={{fontSize:13,color:'#94a3b8'}}>Positions: <span style={{fontFamily:'monospace',color:'#f59e0b'}}>{r.cys.positions.join(', ')||'None'}</span></div>
                  </div>
                  {r.cys.closePairs.length>0?(
                    <>
                      <div style={{fontSize:12,color:'#f87171',marginBottom:8,fontWeight:600}}>⚠ Close Cys pairs (≤6 aa — aggregation / misfolding risk):</div>
                      {r.cys.closePairs.map((p,i)=>(
                        <div key={i} style={{background:'#1c0a0a',border:'1px solid #7f1d1d',borderRadius:6,padding:'8px 12px',marginBottom:6,fontSize:12}}>
                          <span style={{color:'#f87171',fontFamily:'monospace',fontWeight:600}}>Cys{p.pos1} — Cys{p.pos2}</span>
                          <span style={{color:'#78716c',marginLeft:8}}>{p.dist} aa apart</span>
                        </div>
                      ))}
                      <div style={{marginTop:10,fontSize:12,color:'#78716c',lineHeight:1.8,background:'#1c1917',padding:'10px 12px',borderRadius:6,borderLeft:'3px solid #92400e'}}>
                        <strong style={{color:'#fde68a'}}>Recommendations ({orgCfg.label}):</strong><br/>
                        {org==='ecoli'
                          ?<>• Add 1–5 mM DTT or 10 mM 2-ME to all lysis and purification buffers<br/>• Consider SHuffle T7 Express (NEB) if aggregation persists<br/>• Verify any additional cofactor requirements for your protein</>
                          :org==='ylip'
                          ?<>• Y. lipolytica secretory pathway supports disulfide bond formation<br/>• Ensure correct signal peptide for secretion if disulfides are required<br/>• Verify cofactor and media requirements specific to your protein</>
                          :null}
                      </div>
                    </>
                  ):<div style={{color:'#22c55e',fontSize:13}}>✓ No closely spaced Cys pairs. Standard conditions sufficient.</div>}
                </div>
              )}

              {tab==='map'&&(
                <div>
                  <div style={{fontSize:11,color:'#64748b',marginBottom:8,display:'flex',gap:12,flexWrap:'wrap',alignItems:'center'}}>
                    <span>Wrapping codon view. Hover highlighted codons for details.</span>
                    <span style={{background:'#78350f',color:'#fde68a',padding:'1px 6px',borderRadius:2,fontSize:10}}>■ Rare (&lt;10)</span>
                    <span style={{background:'#7f1d1d',color:'#fca5a5',padding:'1px 6px',borderRadius:2,fontSize:10}}>■ Very rare (&lt;5)</span>
                  </div>
                  <CodonMap seq={r.originalSeq} replacements={r.replacements}/>
                </div>
              )}

              {tab==='aa'&&(()=>{
                const maxN=Math.max(...sortedAA.map(a=>r.aaCounts[a.o]||0),1);
                const thStyle=(active)=>({padding:'6px 10px',fontWeight:active?700:600,fontSize:10,letterSpacing:'0.04em',cursor:'pointer',userSelect:'none',color:active?'#fff':'#475569',whiteSpace:'nowrap',border:'none',background:active?'rgba(14,165,233,0.18)':'transparent',borderBottom:active?'2px solid #0ea5e9':'1px solid #1e293b',textAlign:'left',transition:'background 0.15s'});
                return(
                  <div style={{overflowX:'auto'}}>
                    <div style={{fontSize:11,color:'#475569',marginBottom:8}}>Click <strong style={{color:'#94a3b8'}}>Amino Acid</strong> or <strong style={{color:'#94a3b8'}}>%</strong> to sort.</div>
                    <table style={{width:'100%',borderCollapse:'collapse',fontSize:11}}>
                      <thead><tr>
                        <th style={thStyle(!!nameDir)} onClick={()=>toggleSort('name')}>Amino Acid<SortArrow dir={nameDir}/></th>
                        <th style={{...thStyle(false),cursor:'default'}}>3-Letter</th>
                        <th style={{...thStyle(false),cursor:'default'}}>1-Letter</th>
                        <th style={{...thStyle(false),cursor:'default',textAlign:'right'}}>Count</th>
                        <th style={{...thStyle(false),cursor:'default',width:'30%'}}>Frequency</th>
                        <th style={thStyle(!!countDir)} onClick={()=>toggleSort('count')}>%<SortArrow dir={countDir}/></th>
                      </tr></thead>
                      <tbody>{sortedAA.map(aa=>{
                        const n=r.aaCounts[aa.o]||0,pct=s.proteinLength?(n/s.proteinLength*100).toFixed(1):'0.0';
                        const col=AA_COLORS[aa.o]||'#3b82f6';
                        return(
                          <tr key={aa.o} style={{borderBottom:'1px solid #0f172a',opacity:n===0?0.3:1}}>
                            <td style={{padding:'5px 8px',color:'#94a3b8'}}>{aa.f}</td>
                            <td style={{padding:'5px 8px',color:'#64748b',fontFamily:'monospace'}}>{aa.t}</td>
                            <td style={{padding:'5px 8px',fontFamily:'monospace',fontWeight:700,color:col,fontSize:13}}>{aa.o}</td>
                            <td style={{padding:'5px 8px',color:'#e2e8f0',textAlign:'right'}}>{n}</td>
                            <td style={{padding:'5px 8px'}}>
                              <div style={{background:'#1e293b',borderRadius:3,height:10,overflow:'hidden'}}>
                                <div style={{width:`${(n/maxN)*100}%`,background:col,height:'100%',borderRadius:3,transition:'width 0.4s'}}/>
                              </div>
                            </td>
                            <td style={{padding:'5px 8px',color:'#64748b',textAlign:'right',whiteSpace:'nowrap'}}>{pct}%</td>
                          </tr>);
                      })}</tbody>
                    </table>
                  </div>
                );
              })()}

              {tab==='opt'&&(
                <div>
                  <div style={{display:'flex',gap:6,marginBottom:12,flexWrap:'wrap'}}>
                    <CopyBtn text={r.optimizedSeq} label="Copy optimized seq"/>
                    <CopyBtn text={r.protein.replace(/\*$/,'')} label="Copy protein"/>
                  </div>
                  <OptimizedDisplay originalSeq={r.originalSeq} optimizedSeq={r.optimizedSeq}/>
                  <div style={{marginTop:12,fontSize:11,color:'#475569',marginBottom:4}}>Protein sequence:</div>
                  <div style={{fontFamily:'monospace',fontSize:11,color:'#64748b',background:'#0a0f1a',borderRadius:6,padding:'8px 12px',border:'1px solid #1e293b',lineHeight:1.9,wordBreak:'break-all'}}>{r.protein.replace(/\*$/,'')}</div>
                </div>
              )}

            </div>
          </div>

          {/* Export */}
          <div style={{background:'#0f172a',border:'1px solid #1e293b',borderRadius:8,padding:'12px 16px',display:'flex',gap:10,alignItems:'center',flexWrap:'wrap'}}>
            <span style={{fontSize:12,color:'#475569',flex:1}}>Export: <strong style={{color:'#94a3b8'}}>{cur.name}</strong> <span style={{color:'#334155'}}>({orgCfg.label})</span></span>
            <button onClick={()=>dl(toMD(cur.name,r,orgCfg.label),`${cur.name.replace(/[\s/\\]/g,'_')}_report.md`,'text/markdown')}
              style={{padding:'7px 14px',fontSize:12,background:'#1e293b',color:'#94a3b8',border:'1px solid #334155',borderRadius:5,cursor:'pointer'}}>↓ Markdown (.md)</button>
            <button onClick={()=>dl(toWordHTML(cur.name,r,orgCfg.label),`${cur.name.replace(/[\s/\\]/g,'_')}_report.doc`,'application/msword')}
              style={{padding:'7px 14px',fontSize:12,background:'#1e3a5f',color:'#93c5fd',border:'1px solid #1d4ed8',borderRadius:5,cursor:'pointer'}}>↓ Word (.doc)</button>
          </div>
        </>}

        <div style={{textAlign:'center',fontSize:11,color:'#1e293b',paddingTop:4}}>{orgCfg.footer}</div>
      </div>
    </div>
  );
}
