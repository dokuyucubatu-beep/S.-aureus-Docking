# S. aureus DHFR–trimetoprim docking çalışması

Bu repo AChE örneğindeki tekrarlanabilirlik yaklaşımını S. aureus DHFR'ye
uygular. Hedef yapı **2W9H**, kristal ligand **trimetoprim (TOP)**.

- **Lys158 silinmedi:** eksik sekiz ağır atom ve terminal OXT modellendi.
  Gözlenmiş ağır atomlar sabit tutularak eklenen kısım minimize edildi.
- **Exhaustiveness 32:** seed 2026–2035, koşu başına CPU 2, en fazla 20 poz.
- Başlangıç ligandı kristalden alınmadı; SMILES'ten ETKDGv3/MMFF94s ile üretildi.
- Nötr trimetoprim kullanıldı; alternatif protonasyon durumları denenmedi.
- Kristal referansın A alternatif konumu kullanıldı; sular ve EDO çıkarıldı.
- RMSD reseptörün koordinat sisteminde, ligand üst üste bindirilmeden hesaplanır.
- Birinci pozun başarısı ile herhangi bir pozun başarılı olması ayrı raporlanır.

Kurulum ve çalıştırma için ana [README](../README.md) dosyasına bakın.
Hazırlama ayrıntıları [denetim kaydında](../data/prepared/preparation_audit.json),
hesaplama çıktıları `experiments/independent/` altında bulunur.

Bu hesaplama bakteriye karşı etkinlik, MIC veya klinik kullanım kanıtı değildir.
Tek hedef ve tek ligand için poz geri kazanımı ile hesaplamanın tekrarlanabilirliği
incelenir. Modellenen Lys158 koordinatları deneysel gözlem olarak yorumlanmamalıdır.
