from pathlib import Path
import base64,xml.etree.ElementTree as ET, re, json

class Parser:
 def __init__(self,t):self.t=t;self.p=0;self.nodes=[None]
 def number(self):
  c=self.t[self.p];self.p+=1
  if ord(c)<91:return ord(c)-33
  val=0;pow=1
  while ord(c)<107:
   assert 91<=ord(c)<=106,(self.p,c)
   val+=(ord(c)-91)*pow;pow*=16
   c=self.t[self.p];self.p+=1
  assert 107<=ord(c)<=122,(self.p,c)
  return val+(ord(c)-107)*pow
 def node(self):
  c=self.t[self.p];self.p+=1
  if c=='F':
   n=self.number()
   assert 0<n<len(self.nodes),(self.p,n,len(self.nodes))
   return self.nodes[n]
  ob={'type':c,'index':len(self.nodes)};self.nodes.append(ob)
  if c in ('I','%','Q'):
   n=self.number();ob['value']=self.t[self.p:self.p+n];self.p+=n
   if c in ('I','%'):
    assert self.t[self.p]=='G',(self.p,self.t[self.p:self.p+20],ob)
    self.p+=1
   if c!='%':ob['attr']=self.node()
  elif c=='6':
   n=self.number()-1;ob['children']=[self.node() for _ in range(n)]
  elif c in ('-','/'):
   ob['children']=[self.node(),self.node()]
  else:raise ValueError(('unexpected',self.p,c,self.t[max(0,self.p-50):self.p+80]))
  return ob

def value(node):
 if node['type'] in ('I','%','Q'):return node['value'].replace('~',' ')
 return render(node)
def seq(node):return node.get('children',[]) if node['type']=='6' else [node]
def attrmap(args):
 out={}
 for arg in args:
  if arg['type']=='/':out[value(arg['children'][0])]=value(arg['children'][1])
 return out
entities={'&coloneq;':':=','&uminus0;':'-','&srarr;':'->','&InvisibleTimes;':'*','&ApplyFunction;':'','&pi;':'Pi','&le;':'<=','&ge;':'>=','&ne;':'<>','&lt;':'<','&gt;':'>','&equals;':'=','&minus;':'-','&plus;':'+','&comma;':',','&period;':'.','&colon;':':','&semi;':';','&ast;':'*','&sol;':'/','&Hat;':'^','&Lambda;':'Lambda','&Omega;':'Omega','&Phi;':'Phi','&psi;':'psi','&alpha;':'alpha','&sdot;':'*','&#923;':'Lambda','&#937;':'Omega','&#8788;':':=','&#8722;':'-','&Assign;':':=','&ImaginaryI;':'I','&map;':'->'}
def textclean(x):
 x=x.replace('~',' ')
 for a,b in entities.items():x=x.replace(a,b)
 return x
unknown=set()
def render(n):
 typ=n['type']
 if typ in ('I','%','Q'):return textclean(n['value'])
 if typ=='6':return ''.join(render(c) for c in n['children'])
 if typ=='/':return '' # typesetting property, never visible expression equality
 assert typ=='-',typ
 name=value(n['children'][0]);args=seq(n['children'][1]);atts=attrmap(args)
 visible=[a for a in args if a['type']!='/']
 if name in ('mi','mn','mo','mtext'):return render(visible[0]) if visible else ''
 if name=='ms':return '"'+render(visible[0]).replace('"','\\"')+'"'
 if name=='mrow':return ''.join(render(a) for a in visible)
 if name=='mspace':return '\n' if atts.get('linebreak')=='newline' else ''
 if name=='mfenced':return atts.get('open','(')+','.join(render(a) for a in visible)+atts.get('close',')')
 if name=='msup':return '('+render(visible[0])+')^('+render(visible[1])+')'
 if name=='msub':return render(visible[0])+'['+render(visible[1])+']'
 if name=='msubsup':return '('+render(visible[0])+'['+render(visible[1])+'])^('+render(visible[2])+')'
 if name=='mfrac':return '('+render(visible[0])+')/('+render(visible[1])+')'
 if name=='msqrt':return 'sqrt('+''.join(render(a) for a in visible)+')'
 if name=='mroot':return '('+render(visible[0])+')^(1/('+render(visible[1])+'))'
 if name in ('mstyle','mpadded','menclose','semantics'):return ''.join(render(a) for a in visible)
 unknown.add(name)
 return '[['+name+':'+','.join(render(a) for a in visible)+']]'

def extract(path,out):
 root=ET.parse(path).getroot()
 result=[]
 for gi, group in enumerate(root.iter('Group')):
  for kind in ['Input','Output']:
   block=group.find(kind)
   if block is None:continue
   for ei,element in enumerate(block.iter('Equation')):
    raw=(element.attrib.get('display') if kind=='Output' else None) or element.text or ''
    if not raw.strip():continue
    t=base64.b64decode(raw).decode('utf-8')
    p=Parser(t)
    try:
     node=p.node();assert not t[p.p:].strip(),(p.p,t[p.p:p.p+100])
     s=render(node)
     result.append((gi,kind,ei,s))
    except Exception as exc:
     (out/(path.stem+f'_group{gi}_{kind}_{ei}_fail.txt')).write_text(t)
     if kind=='Input': raise RuntimeError((path.name,gi,kind,ei)) from exc
     result.append((gi,kind,ei,'[Non-typeset Maple output omitted; raw saved separately]'))
  # ordinary text outputs contain useful cached export messages
  output=group.find('Output')
  if output is not None:
   for element in output.iter('Text-field'):
    if not any(ch.tag=='Equation' for ch in element.iter()):
     tx=''.join(element.itertext())
     if tx.strip():result.append((gi,'OutputText',0,tx))
 target=out/(path.stem+'.decoded.txt')
 target.write_text('\n\n'.join(f'===== Group {gi} {kind} Equation {ei} =====\n{s}' for gi,kind,ei,s in result))
 return target

if __name__=='__main__':
 import sys
 root=Path(__file__).resolve().parents[2]
 out=root/'code_review/extracted';out.mkdir(parents=True,exist_ok=True)
 paths=[Path(x) for x in sys.argv[1:]] or list((root.parent/'maple').glob('*.mw'))
 for path in paths:print(extract(path,out))
 print('Unknown rendering names:',unknown)
