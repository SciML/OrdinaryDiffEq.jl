struct VTSRKTableau{T,T2}
  d0::T
  d7::T
  n2::T
  n3::T
  n6::T
  n8::T
  q20::T
  q21::T
  q30::T
  q32::T
  q41::T
  q43::T
  q51::T
  q54::T
  q61::T
  q65::T
  q70::T
  q71::T
  q76::T
  q80::T
  q81::T
  q82::T
  q87::T
  r::T
  c2::T2
  c3::T2
  c4::T2
  c5::T2
  c6::T2
  c7::T2
  c8::T2
end

function VTSRKTableau(T,T2)
  d0 = convert(T, 1.0)
  d7 = convert(T, 0.003674184820260)
  n2 = convert(T, 0.179502832154858)
  n3 = convert(T, 0.073789956884809)
  n6 = convert(T, 0.017607159013167)
  n8 = convert(T, 0.729100051947166)
  q20 = convert(T, 0.085330772947643)
  q21 = convert(T, 0.914669227052357)
  q30 = convert(T, 0.058121281984411)
  q32 = convert(T, 0.941878718015589)
  q41 = convert(T, 0.036365639242841)
  q43 = convert(T, 0.802870131352638)
  q51 = convert(T, 0.491214340660555)
  q54 = convert(T, 0.508785659339445)
  q61 = convert(T, 0.566135231631241)
  q65 = convert(T, 0.433864768368758)
  q70 = convert(T, 0.020705281786630)
  q71 = convert(T, 0.091646079651566)
  q76 = convert(T, 0.883974453741544)
  q80 = convert(T, 0.008506650138784)
  q81 = convert(T, 0.110261531523242)
  q82 = convert(T, 0.030113037742445)
  q87 = convert(T, 0.851118780595529)
  r = convert(T, 3.5794)
  c2 = convert(T2, 1/7)
  c3 = convert(T2, 2/7)
  c4 = convert(T2, 3/7)
  c5 = convert(T2, 4/7)
  c6 = convert(T2, 5/7)
  c7 = convert(T2, 6/7)
  c8 = convert(T2, 1.0)
  VTSRKTableau(d0,d7,n2,n3,n6,n8,q20,q21,q30,q32,q41,q43,q51,q54,q61,q65,q70,q71,q76,q80,q81,q82,q87,r,c2,c3,c4,c5,c6,c7,c8)
end
